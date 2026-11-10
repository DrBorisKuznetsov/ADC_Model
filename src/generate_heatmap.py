import os
import numpy as np
import matplotlib.pyplot as plt
from src.solver import simulate_adc_input, SignalGenerator, CapacitorModel, quantize_samples, reconstruct_voltage
from src.analyzer import calculate_fft_metrics

def generate_heatmap():
    plot_dir = "wiki/assets"
    os.makedirs(plot_dir, exist_ok=True)
    
    f_s = 100e3         # 100 кГц
    N_samples = 256     # Снизим для скорости
    T_acq = 300e-9
    R_ext = 200.0
    R_sw = 50.0
    C_sh = 20e-12
    bits = 16
    V_ref = 5.0
    C_nom = 10e-9
    V_rated = 5.0
    
    cap_type = 'X7R'
    cap_model = CapacitorModel(cap_type=cap_type, c_nom=C_nom, v_rated=V_rated)
    
    # Сетка (Частота vs Амплитуда)
    freqs = np.logspace(2, 4.3, 15) # 100 Гц - ~20 кГц
    amps = np.linspace(0.5, 4.5, 15) # 0.5В - 4.5В (от пика до пика)
    
    enob_grid = np.zeros((len(amps), len(freqs)))
    
    print("Генерация тепловой карты...")
    for i, amp in enumerate(amps):
        for j, freq in enumerate(freqs):
            sig_offset = 2.5
            freq_adjusted = np.round(freq / (f_s/N_samples)) * (f_s/N_samples)
            if freq_adjusted == 0: freq_adjusted = f_s/N_samples
            
            sig_gen = SignalGenerator(form='sine', amplitude=amp/2.0, frequency=freq_adjusted, dc_offset=sig_offset)
            
            sim_res = simulate_adc_input(
                signal_gen=sig_gen, cap_model=cap_model,
                R_ext=R_ext, R_sw=R_sw, C_sh=C_sh, T_acq=T_acq, f_s=f_s, N_samples=N_samples, V_ref=V_ref,
                method='semi_implicit', steps_per_cycle=50
            )
            codes, _ = quantize_samples(sim_res['v_samples'], V_ref=V_ref, bits=bits)
            v_rec = reconstruct_voltage(codes, V_ref=V_ref, bits=bits)
            
            metrics = calculate_fft_metrics(
                v_samples=v_rec, f_s=f_s, bits=bits, v_ref=V_ref, v_ext_samples=sim_res['v_ext_samples']
            )
            enob_grid[i, j] = metrics['enob']

    plt.figure(figsize=(10, 8))
    X, Y = np.meshgrid(freqs, amps)
    loss_grid = 16 - enob_grid
    
    plt.pcolormesh(X, Y, loss_grid, shading='auto', cmap='inferno')
    plt.xscale('log')
    plt.colorbar(label='ENOB Loss (bits)')
    plt.title('Heatmap: ENOB Degradation due to X7R Capacitor Nonlinearity\n(100 kS/s Sample Rate, Ideal 16-bit ADC)')
    plt.xlabel('Input Signal Frequency (Hz)')
    plt.ylabel('Signal Amplitude (V p-p)')
    
    CS = plt.contour(X, Y, loss_grid, levels=[1, 3, 5, 7], colors='white', alpha=0.5)
    plt.clabel(CS, inline=True, fontsize=10, fmt='%1.1f bits')
    
    output_path = os.path.join(plot_dir, "enob_heatmap.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Heatmap saved to {output_path}")

if __name__ == '__main__':
    generate_heatmap()
