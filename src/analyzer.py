import numpy as np
import scipy.fft as fft
from typing import Dict, Any, List, Tuple, Optional
from src.solver import simulate_adc_input, SignalGenerator, CapacitorModel, quantize_samples, reconstruct_voltage

def calculate_fft_metrics(
    v_samples: np.ndarray,
    f_s: float,
    bits: int,
    v_ref: float,
    v_ext_samples: Optional[np.ndarray] = None
) -> Dict[str, Any]:
    """
    Вычисляет спектральные метрики (THD, SINAD, ENOB, SNR) и ошибку установления.
    
    Параметры:
        v_samples: массив аналоговых отсчетов на входе АЦП (В)
        f_s: частота дискретизации (Гц)
        bits: разрядность АЦП
        v_ref: опорное напряжение (В)
        v_ext_samples: (опционально) массив внешних напряжений фильтра в моменты выборки
        
    Возвращает:
        Словарь с метриками:
            'thd_db': THD в дБ
            'sinad_db': SINAD в дБ
            'enob': эффективная разрядность (биты)
            'enob_loss': потери эффективной разрядности (биты)
            'snr_db': SNR в дБ
            'settling_error_lsb': максимальная ошибка установления в LSB
            'frequencies': массив частот для спектрограммы
            'power_spec_db': спектр мощности в дБ (нормированный к пику)
            'harmonics_idx': индексы найденных гармоник в спектре
    """
    N = len(v_samples)
    
    # 1. Применяем оконную функцию Блэкмана-Харриса (4-термовую)
    n_arr = np.arange(N)
    window = 0.35875 - 0.48829 * np.cos(2 * np.pi * n_arr / (N - 1)) + 0.14128 * np.cos(4 * np.pi * n_arr / (N - 1)) - 0.01168 * np.cos(6 * np.pi * n_arr / (N - 1))
    
    # Расчет коэффициентов когерентного усиления и шумовой полосы окна
    cg = np.mean(window)
    npg = np.mean(window ** 2)
    # Коррекция соотношения сигнал/шум для оконного БПФ
    window_correction = cg ** 2 / npg
    
    # Remove DC with the same window weights used by the FFT. An ordinary
    # arithmetic mean reintroduces low-frequency energy after windowing when
    # the record contains a non-integer number of signal periods.
    weighted_mean = np.sum(v_samples * window) / np.sum(window)
    v_windowed = (v_samples - weighted_mean) * window
    
    # БПФ (для вещественного сигнала)
    fft_val = fft.rfft(v_windowed)
    mag = np.abs(fft_val)
    power = mag ** 2
    
    # Защита от деления на ноль в пустом спектре
    if np.sum(power) == 0:
        power += 1e-20
    
    # Нормируем частоты
    freqs = fft.rfftfreq(N, 1.0 / f_s)
    
    # Находим индекс фундаментальной частоты, исключая только DC-бин.
    # Низкочастотный сигнал может попадать в первые FFT-бины; их нельзя
    # выбрасывать как шум около нуля, иначе SINAD/ENOB будут завышены.
    analysis_start = 1
    if len(power) > analysis_start:
        fund_idx = np.argmax(power[analysis_start:]) + analysis_start
    else:
        fund_idx = 0
        
    # Ширина пика фундаментальной частоты для окна Блэкмана-Харриса (берем 11 бинов: +-5)
    half_width = 5
    fund_bins = np.arange(max(analysis_start, fund_idx - half_width), min(len(power), fund_idx + half_width + 1))
    
    # Мощность основного сигнала в спектре
    P_fund = np.sum(power[fund_bins])
    
    # Нахождение гармоник (с учетом алиасинга/наложения спектров)
    num_harmonics = 10
    P_harmonics_list = []
    harmonics_idx_list = []
    
    for k in range(2, num_harmonics + 1):
        h_idx_unaliased = k * fund_idx
        N_fft = len(power)
        N_orig = 2 * (N_fft - 1)
        h_idx = h_idx_unaliased % N_orig
        if h_idx >= N_fft:
            h_idx = N_orig - h_idx
            
        h_bins = np.arange(max(0, h_idx - half_width), min(len(power), h_idx + half_width + 1))
        # Исключаем наложение на основной сигнал и DC-бин.
        h_bins = np.array([b for b in h_bins if b not in fund_bins and b >= analysis_start])
        
        if len(h_bins) > 0:
            P_h = np.sum(power[h_bins])
            P_harmonics_list.append(P_h)
            harmonics_idx_list.append(h_idx)
        else:
            P_harmonics_list.append(0.0)
            
    P_harm_total = np.sum(P_harmonics_list)
    
    # Общая мощность анализируемого спектра (исключая только DC-бин).
    P_total_band = np.sum(power[analysis_start:])
    
    # Мощность шума и искажений (все кроме основного сигнала в анализируемой полосе)
    P_nd = max(1e-25, P_total_band - P_fund)
    
    # Мощность чистого шума (без гармоник)
    P_noise_only = max(1e-25, P_nd - P_harm_total)
    
    # Применяем коррекцию окна к мощностям шума и гармоник перед расчетом логарифмов
    P_nd_corrected = P_nd / window_correction
    P_harm_total_corrected = P_harm_total / window_correction
    P_noise_only_corrected = P_noise_only / window_correction
    
    # Вычисление THD
    if P_fund > 0 and P_harm_total_corrected > 0:
        thd_db = 10 * np.log10(P_harm_total_corrected / P_fund)
    else:
        thd_db = -120.0
        
    # SINAD
    sinad_corrected = P_fund / P_nd_corrected
    sinad_db = 10 * np.log10(sinad_corrected)
    
    # SNR
    snr_db = 10 * np.log10(P_fund / P_noise_only_corrected)
    
    # Расчет ENOB. Физически эффективная разрядность не должна превышать
    # выбранную номинальную разрядность АЦП.
    enob_raw = (sinad_db - 1.76) / 6.02
    enob = min(float(bits), max(0.0, enob_raw))
    enob_loss = max(0.0, bits - enob)
    
    # Нормализованный спектр мощности в дБ для графиков (0 дБ = пик основного сигнала)
    max_power = max(1e-25, np.max(power))
    power_db = 10 * np.log10(np.maximum(power / max_power, 1e-14))
    power_db = np.clip(power_db, -140.0, 0.0)
    
    # Расчет максимальной ошибки установления в LSB (исключая переходный процесс запуска)
    settling_error_lsb = 0.0
    if v_ext_samples is not None:
        lsb = v_ref / ((2 ** bits) - 1)
        # Исключаем первые 10 отсчетов для завершения переходного процесса запуска
        discard = min(len(v_samples) // 10, 10)
        errors = np.abs(v_ext_samples[discard:] - v_samples[discard:])
        settling_error_lsb = np.max(errors) / lsb
        
    return {
        'thd_db': float(thd_db),
        'sinad_db': float(sinad_db),
        'enob_raw': float(enob_raw),
        'enob': float(enob),
        'enob_loss': float(enob_loss),
        'snr_db': float(snr_db),
        'settling_error_lsb': float(settling_error_lsb),
        'frequencies': freqs,
        'power_spec_db': power_db,
        'harmonics_idx': harmonics_idx_list
    }

def assess_spectral_record(
    signal_frequency: float,
    f_s: float,
    n_samples: int,
    signal_form: str = 'sine'
) -> Dict[str, Any]:
    """Assess whether the current record supports a sine-based ENOB estimate."""
    if signal_frequency <= 0 or f_s <= 0 or n_samples <= 0:
        raise ValueError("Signal frequency, sampling rate and sample count must be positive")

    alias_frequency = abs(((signal_frequency + 0.5 * f_s) % f_s) - 0.5 * f_s)
    observed_cycles = alias_frequency * n_samples / f_s
    bin_width = f_s / n_samples
    warnings: List[str] = []

    if signal_form.lower() != 'sine':
        warnings.append("ENOB/SINAD assumes a sine input; waveform harmonics are counted as distortion")
    if alias_frequency <= max(1e-12, 1e-12 * f_s):
        warnings.append("Input aliases to DC at the selected sampling rate")
    elif observed_cycles < 4.0:
        warnings.append(
            f"Only {observed_cycles:.2f} aliased signal cycles are present; use at least 4 cycles"
        )
    if 0.5 * f_s - alias_frequency < bin_width:
        warnings.append("Input is within one FFT bin of Nyquist; the result is phase-sensitive")

    return {
        'valid_for_enob': not warnings,
        'alias_frequency_hz': float(alias_frequency),
        'observed_cycles': float(observed_cycles),
        'warnings': warnings
    }

def run_parametric_sweep(
    base_params: Dict[str, Any],
    sweep_axes: Dict[str, List[Any]]
) -> Dict[str, Any]:
    """
    Многопараметрический анализ (N-Dimensional Space Sweep).
    
    Параметры:
        base_params: базовые параметры симулятора (R_ext, R_sw, C_ext_nom, C_sh, etc.)
        sweep_axes: оси сканирования, например:
                    {
                        'R_ext': [10.0, 100.0, 1000.0],
                        'sig_frequency': [1000.0, 10000.0]
                    }
                    
    Возвращает:
        Словарь с осями сканирования и результатами (гиперматрицами ENOB, THD, clipping)
    """
    import itertools
    
    # Извлекаем имена и списки значений сканируемых параметров
    axis_names = list(sweep_axes.keys())
    axis_values = [sweep_axes[name] for name in axis_names]
    
    shape = [len(vals) for vals in axis_values]
    
    # Создаем результирующие многомерные массивы
    sweep_enob = np.zeros(shape)
    sweep_enob_loss = np.zeros(shape)
    sweep_thd = np.zeros(shape)
    sweep_clipping = np.zeros(shape, dtype=bool)
    
    direct_param_names = {
        'R_ext', 'R_sw', 'C_sh', 'T_acq', 'f_s', 'N_samples', 'V_ref', 'method', 'bits',
        'sig_form', 'sig_amplitude', 'sig_frequency', 'sig_dc_offset',
        'cap_c_nom', 'cap_v_rated', 'cap_cap_type', 'cap_custom_drop', 'cap_v50', 'cap_n'
    }
    legacy_param_aliases = {
        'form': 'sig_form',
        'amplitude': 'sig_amplitude',
        'frequency': 'sig_frequency',
        'dc_offset': 'sig_dc_offset',
        'c_nom': 'cap_c_nom',
        'v_rated': 'cap_v_rated',
        'cap_type': 'cap_cap_type',
        'custom_drop': 'cap_custom_drop',
        'v50': 'cap_v50',
        'n': 'cap_n',
    }
    
    index_ranges = [range(len(vals)) for vals in axis_values]
    for multi_idx in itertools.product(*index_ranges):
        combo = tuple(axis_values[i][multi_idx[i]] for i in range(len(axis_names)))
        run_params = base_params.copy()
        
        # Обновляем параметры для запуска
        for i, name in enumerate(axis_names):
            target_name = legacy_param_aliases.get(name, name)
            if target_name not in direct_param_names:
                raise ValueError(f"Unsupported sweep axis: {name}")
            run_params[target_name] = combo[i]
                
        # Создаем объекты генератора и модели конденсатора
        sig_gen = SignalGenerator(
            form=run_params.get('sig_form', 'sine'),
            amplitude=run_params.get('sig_amplitude', 4.0),
            frequency=run_params.get('sig_frequency', 1000.0),
            dc_offset=run_params.get('sig_dc_offset', 5.0)
        )
        
        cap_model = CapacitorModel(
            cap_type=run_params.get('cap_cap_type', 'X7R'),
            c_nom=run_params.get('cap_c_nom', 10e-9),
            v_rated=run_params.get('cap_v_rated', 10.0),
            custom_drop=run_params.get('cap_custom_drop', None),
            v50=run_params.get('cap_v50', None),
            n=run_params.get('cap_n', None)
        )
        
        # Симулируем
        res = simulate_adc_input(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=run_params.get('R_ext', 100.0),
            R_sw=run_params.get('R_sw', 50.0),
            C_sh=run_params.get('C_sh', 20e-12),
            T_acq=run_params.get('T_acq', 1e-7),
            f_s=run_params.get('f_s', 1e5),
            N_samples=run_params.get('N_samples', 128),
            V_ref=run_params.get('V_ref', 10.0),
            method=run_params.get('method', 'semi_implicit'),
            steps_per_cycle=25,
            store_trace=False
        )
        
        # Квантуем
        codes, clip = quantize_samples(
            res['v_samples'], 
            V_ref=run_params.get('V_ref', 10.0), 
            bits=run_params.get('bits', 16)
        )
        v_rec = reconstruct_voltage(codes, V_ref=run_params.get('V_ref', 10.0), bits=run_params.get('bits', 16))
        
        # Анализируем спектр
        metrics = calculate_fft_metrics(
            v_samples=v_rec,
            f_s=run_params.get('f_s', 1e5),
            bits=run_params.get('bits', 16),
            v_ref=run_params.get('V_ref', 10.0)
        )
        
        # Сохраняем результаты
        sweep_enob[multi_idx] = metrics['enob']
        sweep_enob_loss[multi_idx] = metrics['enob_loss']
        sweep_thd[multi_idx] = metrics['thd_db']
        sweep_clipping[multi_idx] = clip
        
    return {
        'axes': sweep_axes,
        'axis_names': axis_names,
        'enob': sweep_enob,
        'enob_loss': sweep_enob_loss,
        'thd': sweep_thd,
        'clipping': sweep_clipping
    }
