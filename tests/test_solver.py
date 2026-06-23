import unittest
import warnings
import numpy as np
from unittest.mock import patch
from src.solver import (
    CAPACITOR_PRESETS,
    SignalGenerator,
    CapacitorModel,
    simulate_adc_input,
    quantize_samples,
    reconstruct_voltage,
)
from src.analyzer import calculate_fft_metrics, run_parametric_sweep, assess_spectral_record

class TestSignalGenerator(unittest.TestCase):
    def test_sine_wave(self):
        gen = SignalGenerator(form='sine', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        # В моменты t = 0, 1/2000, 1/1000
        self.assertAlmostEqual(gen(0.0), 2.5)
        self.assertAlmostEqual(gen(0.25 / 1000.0), 4.5)  # пик синуса
        self.assertAlmostEqual(gen(0.5 / 1000.0), 2.5)
        self.assertAlmostEqual(gen(0.75 / 1000.0), 0.5)  # впадина синуса

    def test_square_wave(self):
        gen = SignalGenerator(form='square', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        self.assertAlmostEqual(gen(0.1 / 1000.0), 4.5)
        self.assertAlmostEqual(gen(0.6 / 1000.0), 0.5)

    def test_triangle_wave(self):
        gen = SignalGenerator(form='triangle', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        self.assertAlmostEqual(gen(0.0), 2.5)
        self.assertAlmostEqual(gen(0.25 / 1000.0), 4.5)
        self.assertAlmostEqual(gen(0.5 / 1000.0), 2.5)
        self.assertAlmostEqual(gen(0.75 / 1000.0), 0.5)

class TestCapacitorModel(unittest.TestCase):
    def test_shared_capacitor_catalog_units(self):
        preset = CAPACITOR_PRESETS['murata_grm155_10u_10v']
        self.assertAlmostEqual(preset['c_nom'], 10e-6)
        self.assertEqual(preset['cap_type'], 'PRESET_FIT')

    def test_c0g_linear(self):
        cap = CapacitorModel(cap_type='C0G', c_nom=10e-9, v_rated=10.0)
        self.assertAlmostEqual(cap.get_C(0.0), 10e-9)
        self.assertAlmostEqual(cap.get_C(5.0), 10e-9)
        self.assertAlmostEqual(cap.get_C(-5.0), 10e-9)
        self.assertAlmostEqual(cap.get_C(20.0), 10e-9) # За границей номинала

    def test_x7r_nonlinear(self):
        cap = CapacitorModel(cap_type='X7R', c_nom=10e-9, v_rated=10.0)
        # В нуле емкость равна номинальной
        self.assertAlmostEqual(cap.get_C(0.0), 10e-9)
        # При номинальном напряжении (10В) емкость должна упасть до 50%
        self.assertAlmostEqual(cap.get_C(10.0), 5e-9)
        # Характеристика симметрична
        self.assertAlmostEqual(cap.get_C(-10.0), 5e-9)
        # Проверка отсечки/ограничения на больших напряжениях (за пределами LUT)
        self.assertAlmostEqual(cap.get_C(30.0), cap.get_C(20.0))

    def test_custom_drop(self):
        # Падение емкости на 40% при номинальном напряжении
        cap = CapacitorModel(cap_type='CUSTOM', c_nom=10e-9, v_rated=10.0, custom_drop=0.4)
        self.assertAlmostEqual(cap.get_C(0.0), 10e-9)
        # При 10В емкость должна быть 6e-9 (падение на 40%)
        self.assertAlmostEqual(cap.get_C(10.0), 6e-9)

    def test_preset_fit(self):
        # Аналитический фитинг Murata GRM155 (10uF, 10V, X5R): v50=3.0, n=1.6
        cap = CapacitorModel(cap_type='PRESET_FIT', c_nom=10e-6, v_rated=10.0, v50=3.0, n=1.6)
        self.assertAlmostEqual(cap.get_C(0.0), 10e-6)
        # При 3.0В емкость должна упасть ровно в 2 раза (до 5uF)
        self.assertAlmostEqual(cap.get_C(3.0), 5e-6)
        self.assertAlmostEqual(cap.get_C(-3.0), 5e-6)
        # При 10В емкость должна быть 10e-6 / (1 + (10/3)**1.6)
        expected_c_10v = 10e-6 / (1.0 + (10.0 / 3.0) ** 1.6)
        self.assertAlmostEqual(cap.get_C(10.0), expected_c_10v)

class TestQuantizer(unittest.TestCase):
    def test_quantize(self):
        v_samples = np.array([-1.0, 0.0, 2.5, 5.0, 6.0])
        codes, clip = quantize_samples(v_samples, V_ref=5.0, bits=8)
        
        # -1.0 должно клиппироваться в 0 -> код 0
        # 6.0 должно клиппироваться в 5.0 -> код 255
        # 2.5 должно быть ровно посередине -> код 128
        self.assertEqual(codes[0], 0)
        self.assertEqual(codes[1], 0)
        self.assertEqual(codes[2], 128)
        self.assertEqual(codes[3], 255)
        self.assertEqual(codes[4], 255)
        
        # Клиппинг должен быть обнаружен
        self.assertTrue(clip)

    def test_reconstruct(self):
        codes = np.array([0, 128, 255])
        v_rec = reconstruct_voltage(codes, V_ref=5.0, bits=8)
        self.assertAlmostEqual(v_rec[0], 0.0)
        self.assertAlmostEqual(v_rec[1], 128.0 * (5.0 / 255.0))
        self.assertAlmostEqual(v_rec[2], 5.0)

class TestSolverIntegration(unittest.TestCase):
    def test_simulation_run(self):
        sig_gen = SignalGenerator(form='sine', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        cap_model = CapacitorModel(cap_type='X7R', c_nom=10e-9, v_rated=5.0)
        
        # Запуск быстрого метода
        res_si = simulate_adc_input(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=100.0,
            R_sw=50.0,
            C_sh=20e-12,
            T_acq=200e-9,
            f_s=100e3,
            N_samples=64,
            V_ref=5.0,
            method='semi_implicit',
            steps_per_cycle=50
        )
        
        # Проверяем размеры массивов
        self.assertEqual(len(res_si['t_samples']), 64)
        self.assertEqual(len(res_si['v_samples']), 64)
        # Суммарное количество шагов: 64 периода по 50 шагов
        self.assertEqual(len(res_si['t']), 64 * 50)
        
        # Значения напряжений должны лежать в физических пределах
        self.assertTrue(np.all(res_si['v_ext'] >= 0.0) and np.all(res_si['v_ext'] <= 5.0))
        self.assertTrue(np.all(res_si['v_sh'] >= 0.0) and np.all(res_si['v_sh'] <= 5.0))

    def test_rejects_acquisition_time_longer_than_sample_period(self):
        sig_gen = SignalGenerator(form='sine', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        cap_model = CapacitorModel(cap_type='X7R', c_nom=10e-9, v_rated=5.0)

        with self.assertRaises(ValueError):
            simulate_adc_input(
                signal_gen=sig_gen,
                cap_model=cap_model,
                R_ext=100.0,
                R_sw=50.0,
                C_sh=20e-12,
                T_acq=20e-6,
                f_s=100e3,
                N_samples=64,
                V_ref=5.0,
                method='semi_implicit',
                steps_per_cycle=50
            )

    def test_simulation_without_trace_keeps_samples(self):
        sig_gen = SignalGenerator(form='sine', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        cap_model = CapacitorModel(cap_type='X7R', c_nom=10e-9, v_rated=5.0)
        common = dict(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=100.0,
            R_sw=50.0,
            C_sh=20e-12,
            T_acq=200e-9,
            f_s=100e3,
            N_samples=64,
            V_ref=5.0,
            method='semi_implicit',
            steps_per_cycle=50,
        )

        with_trace = simulate_adc_input(**common)
        without_trace = simulate_adc_input(**common, store_trace=False)

        self.assertEqual(len(without_trace['t']), 0)
        np.testing.assert_allclose(without_trace['v_samples'], with_trace['v_samples'])
        np.testing.assert_allclose(without_trace['v_ext_samples'], with_trace['v_ext_samples'])

    def test_trace_until_limits_only_plot_trace(self):
        sig_gen = SignalGenerator(form='sine', amplitude=2.0, frequency=1000.0, dc_offset=2.5)
        cap_model = CapacitorModel(cap_type='X7R', c_nom=10e-9, v_rated=5.0)
        common = dict(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=100.0,
            R_sw=50.0,
            C_sh=20e-12,
            T_acq=200e-9,
            f_s=100e3,
            N_samples=64,
            V_ref=5.0,
            method='semi_implicit',
            steps_per_cycle=50,
        )

        full_trace = simulate_adc_input(**common)
        short_trace = simulate_adc_input(**common, trace_until=10e-6)

        self.assertLess(len(short_trace['t']), len(full_trace['t']))
        np.testing.assert_allclose(short_trace['v_samples'], full_trace['v_samples'])

    def test_fast_solver_tracks_radau_at_short_acquisition_edge(self):
        sig_gen = SignalGenerator(form='sine', amplitude=2.2, frequency=20e3, dc_offset=2.5)
        cap_model = CapacitorModel(cap_type='X7R', c_nom=0.1e-9, v_rated=5.0)
        common = dict(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=2000.0,
            R_sw=50.0,
            C_sh=20e-12,
            T_acq=300e-9,
            f_s=100e3,
            N_samples=32,
            V_ref=5.0,
            store_trace=False,
        )

        fast = simulate_adc_input(**common, method='semi_implicit', steps_per_cycle=25)
        reference = simulate_adc_input(**common, method='radau')
        peak_error = np.max(np.abs(fast['v_samples'][10:] - reference['v_samples'][10:]))

        self.assertLess(peak_error, 5e-3)

class TestAnalyzer(unittest.TestCase):
    def test_calculate_fft_metrics(self):
        fs = 100e3
        fin = 1.01e3
        N = 256
        t = np.arange(N) / fs
        v = 2.5 + 2.0 * np.sin(2 * np.pi * fin * t)
        
        metrics = calculate_fft_metrics(v, fs, bits=16, v_ref=5.0)
        
        self.assertIn('thd_db', metrics)
        self.assertIn('sinad_db', metrics)
        self.assertIn('enob', metrics)
        self.assertIn('enob_loss', metrics)
        self.assertIn('snr_db', metrics)
        self.assertIn('settling_error_lsb', metrics)
        self.assertGreater(metrics['enob'], 10.0)

    def test_low_frequency_fundamental_does_not_overstate_enob(self):
        fs = 100e3
        fin = 101.0
        bits = 16
        v_ref = 5.0
        N = 1024
        t = np.arange(N) / fs
        v = 2.7 + 2.2 * np.sin(2 * np.pi * fin * t)

        codes, _ = quantize_samples(v, V_ref=v_ref, bits=bits)
        v_rec = reconstruct_voltage(codes, V_ref=v_ref, bits=bits)
        metrics = calculate_fft_metrics(v_rec, fs, bits=bits, v_ref=v_ref)

        self.assertLess(metrics['sinad_db'], 120.0)
        self.assertGreater(metrics['enob'], 12.0)
        self.assertLessEqual(metrics['enob'], float(bits))
        self.assertGreater(metrics['enob_loss'], 0.0)

    def test_noncoherent_sine_does_not_create_large_false_enob_loss(self):
        fs = 100e3
        fin = 1010.0
        bits = 16
        v_ref = 5.0
        N = 1024
        t = np.arange(N) / fs
        v = 2.5 + 2.2 * np.sin(2 * np.pi * fin * t)

        codes, _ = quantize_samples(v, V_ref=v_ref, bits=bits)
        v_rec = reconstruct_voltage(codes, V_ref=v_ref, bits=bits)
        metrics = calculate_fft_metrics(v_rec, fs, bits=bits, v_ref=v_ref)

        self.assertGreater(metrics['enob'], 13.0)
        self.assertLess(metrics['enob_loss'], 3.0)

    def test_spectral_record_flags_edge_conditions(self):
        low_frequency = assess_spectral_record(10.0, 100e3, 1024, 'sine')
        aliased_dc = assess_spectral_record(100e3, 100e3, 1024, 'sine')
        nonsine = assess_spectral_record(1e3, 100e3, 1024, 'square')
        valid = assess_spectral_record(1010.0, 100e3, 1024, 'sine')

        self.assertFalse(low_frequency['valid_for_enob'])
        self.assertFalse(aliased_dc['valid_for_enob'])
        self.assertFalse(nonsine['valid_for_enob'])
        self.assertTrue(valid['valid_for_enob'])

    def test_constant_record_produces_finite_spectrum_without_warnings(self):
        samples = np.full(128, 2.5)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            metrics = calculate_fft_metrics(samples, 100e3, bits=16, v_ref=5.0)

        self.assertEqual(caught, [])
        self.assertTrue(np.all(np.isfinite(metrics['power_spec_db'])))

    def test_parametric_sweep_uses_ui_signal_frequency_axis(self):
        seen_frequencies = []

        def fake_simulate(signal_gen, cap_model, **kwargs):
            seen_frequencies.append(signal_gen.frequency)
            return {
                'v_samples': np.array([0.1, 0.2, 0.3, 0.4]),
                'v_ext_samples': np.array([0.1, 0.2, 0.3, 0.4]),
            }

        fake_metrics = {
            'enob': 12.0,
            'enob_loss': 4.0,
            'thd_db': -80.0,
        }
        base_params = {
            'R_ext': 100.0,
            'R_sw': 50.0,
            'C_sh': 20e-12,
            'T_acq': 200e-9,
            'f_s': 100e3,
            'N_samples': 16,
            'V_ref': 5.0,
            'bits': 16,
            'method': 'semi_implicit',
            'sig_form': 'sine',
            'sig_amplitude': 2.0,
            'sig_frequency': 1000.0,
            'sig_dc_offset': 2.5,
            'cap_cap_type': 'X7R',
            'cap_c_nom': 10e-9,
            'cap_v_rated': 5.0,
            'cap_custom_drop': None,
            'cap_v50': None,
            'cap_n': None,
        }

        with patch('src.analyzer.simulate_adc_input', side_effect=fake_simulate), \
             patch('src.analyzer.calculate_fft_metrics', return_value=fake_metrics):
            sweep = run_parametric_sweep(
                base_params,
                {'sig_frequency': [500.0, 1000.0, 5000.0]}
            )

        self.assertEqual(seen_frequencies, [500.0, 1000.0, 5000.0])
        self.assertEqual(sweep['enob_loss'].shape, (3,))

if __name__ == '__main__':
    unittest.main()
