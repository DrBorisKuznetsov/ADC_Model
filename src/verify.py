import os
import numpy as np
import matplotlib.pyplot as plt
from src.solver import simulate_adc_input, SignalGenerator, CapacitorModel, quantize_samples, reconstruct_voltage
from src.analyzer import calculate_fft_metrics

def run_verification():
    # Создаем директорию для графиков, если она не существует
    plot_dir = "wiki/assets"
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Задаем параметры симуляции
    f_s = 100e3         # 100 кГц частота дискретизации
    N_samples = 1024    # Число отсчетов для БПФ
    T_acq = 300e-9      # 300 нс время выборки
    R_ext = 200.0       # 200 Ом сопротивление фильтра
    R_sw = 50.0         # 50 Ом сопротивление ключа выборки
    C_sh = 20e-12       # 20 пФ емкость выборки АЦП
    bits = 16           # 16-битный АЦП
    V_ref = 5.0         # 5В опорное напряжение
    
    # Параметры входного сигнала (Синусоида)
    sig_freq = 1.01e3   # ~1 кГц (чуть смещенная, чтобы не совпадать с бинами БПФ идеально)
    sig_amp = 2.2       # Амплитуда 2.2 В (размах 4.4 В)
    sig_offset = 2.5    # Постоянное смещение 2.5 В (сигнал от 0.3В до 4.7В)
    
    # Параметры конденсатора
    C_nom = 10e-9       # 10 нФ номинальная емкость
    V_rated = 5.0       # 5В номинальное напряжение
    
    # Генератор сигнала
    sig_gen = SignalGenerator(form='sine', amplitude=sig_amp, frequency=sig_freq, dc_offset=sig_offset)
    
    # Типы диэлектриков для тестирования
    dielectrics = ['C0G', 'X7R', 'X5R']
    results = {}
    
    print("=" * 60)
    print("ВЕРИФИКАЦИЯ ФИЗИЧЕСКОЙ КОРРЕКТНОСТИ МОДЕЛИ ВХОДНОГО ТРАКТА АЦП")
    print("=" * 60)
    print(f"Частота сигнала: {sig_freq/1e3:.3f} кГц, Амплитуда: {sig_amp} В, Смещение: {sig_offset} В")
    print(f"Параметры фильтра: R_ext={R_ext} Ом, C_ext_nom={C_nom*1e9} нФ (V_rated={V_rated} В)")
    print(f"Параметры АЦП: R_sw={R_sw} Ом, C_sh={C_sh*1e12} пФ, T_acq={T_acq*1e9} нс, Разрядность: {bits} бит")
    print("-" * 60)
    
    # Создаем большой общий график сравнения характеристик C(V) и спектров
    plt.figure(figsize=(15, 10))
    
    # 1-й подграфик: Зависимость C(V)
    plt.subplot(2, 2, 1)
    v_axis = np.linspace(-10, 10, 200)
    for cap_type in dielectrics:
        model = CapacitorModel(cap_type=cap_type, c_nom=C_nom, v_rated=V_rated)
        plt.plot(v_axis, model.get_C(v_axis) * 1e9, label=cap_type, lw=2)
    plt.title("Характеристики конденсаторов C(V)")
    plt.xlabel("Напряжение на конденсаторе (В)")
    plt.ylabel("Емкость (нФ)")
    plt.grid(True)
    plt.legend()
    
    for cap_type in dielectrics:
        # Модель конденсатора
        cap_model = CapacitorModel(cap_type=cap_type, c_nom=C_nom, v_rated=V_rated)
        
        # Симуляция (используем быстрый полунеявный метод)
        sim_res = simulate_adc_input(
            signal_gen=sig_gen,
            cap_model=cap_model,
            R_ext=R_ext,
            R_sw=R_sw,
            C_sh=C_sh,
            T_acq=T_acq,
            f_s=f_s,
            N_samples=N_samples,
            V_ref=V_ref,
            method='semi_implicit',
            steps_per_cycle=100
        )
        
        # Оцифровка
        codes, clip = quantize_samples(sim_res['v_samples'], V_ref=V_ref, bits=bits)
        v_rec = reconstruct_voltage(codes, V_ref=V_ref, bits=bits)
        
        # Спектральный анализ
        metrics = calculate_fft_metrics(
            v_samples=v_rec,
            f_s=f_s,
            bits=bits,
            v_ref=V_ref,
            v_ext_samples=sim_res['v_ext_samples']
        )
        
        results[cap_type] = {
            'metrics': metrics,
            'sim': sim_res,
            'clip': clip
        }
        
        print(f"Диэлектрик: {cap_type}")
        print(f"  ENOB: {metrics['enob']:.2f} бит (Потеря: {metrics['enob_loss']:.2f} бит)")
        print(f"  THD:  {metrics['thd_db']:.2f} дБ")
        print(f"  SNR:  {metrics['snr_db']:.2f} дБ")
        print(f"  Клиппинг обнаружен: {'ДА' if clip else 'НЕТ'}")
        print("-" * 60)
        
    # Подграфики спектров
    plt.subplot(2, 2, 2)
    # Строим спектры для каждого диэлектрика
    for cap_type in dielectrics:
        metrics = results[cap_type]['metrics']
        plt.plot(metrics['frequencies']/1e3, metrics['power_spec_db'], label=f"{cap_type} (ENOB: {metrics['enob']:.1f})", alpha=0.7)
    plt.title("Сравнение спектров мощности FFT")
    plt.xlabel("Частота (кГц)")
    plt.ylabel("Мощность (дБ)")
    plt.xlim(0, f_s / 2 / 1e3)
    plt.ylim(-140, 5)
    plt.grid(True)
    plt.legend()
    
    # 3-й подграфик: Временная область (увеличим первые несколько циклов)
    plt.subplot(2, 1, 2)
    sim_c0g = results['C0G']['sim']
    sim_x5r = results['X5R']['sim']
    
    # Покажем только первые 5 циклов
    n_cycles_plot = 5
    points_to_plot = int(n_cycles_plot * 100) # 100 шагов на цикл
    
    t_ms = sim_c0g['t'][:points_to_plot] * 1e3
    plt.plot(t_ms, sim_c0g['v_in'][:points_to_plot], 'k--', label="V_in (Вход)", alpha=0.7)
    plt.plot(t_ms, sim_c0g['v_ext'][:points_to_plot], 'g-', label="V_ext (C0G - линейный)", lw=2)
    plt.plot(t_ms, sim_x5r['v_ext'][:points_to_plot], 'r-', label="V_ext (X5R - нелинейный)", lw=2)
    
    # Точки выборок
    t_s_ms = sim_c0g['t_samples'][:n_cycles_plot] * 1e3
    plt.scatter(t_s_ms, results['C0G']['sim']['v_samples'][:n_cycles_plot], color='green', marker='o', s=80, zorder=5, label="Отсчеты C0G")
    plt.scatter(t_s_ms, results['X5R']['sim']['v_samples'][:n_cycles_plot], color='red', marker='x', s=80, zorder=5, label="Отсчеты X5R")
    
    # Отрисовка интервалов выборки (is_closed)
    # Подсветим серым фоном зоны выборки
    for n in range(n_cycles_plot):
        plt.axvspan(n * (1.0/f_s)*1e3, (n * (1.0/f_s) + T_acq)*1e3, color='blue', alpha=0.1, label="T_acq (Выборка)" if n==0 else "")
        
    plt.title("Процессы во временной области (первые 5 периодов дискретизации)")
    plt.xlabel("Время (мс)")
    plt.ylabel("Напряжение (В)")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    
    # Сохраняем график
    output_path = os.path.join(plot_dir, "verification_results.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    
    print(f"График верификации сохранен по пути: {output_path}")
    print("=" * 60)

if __name__ == '__main__':
    run_verification()
