import json
import math
from bisect import bisect_left
from pathlib import Path
from typing import Dict, Tuple, Union, Optional

import numpy as np
import scipy.integrate as integrate

class SignalGenerator:
    """Генератор сигналов: Синус, Меандр (Ступеньки), Треугольник."""
    def __init__(self, form: str, amplitude: float, frequency: float, dc_offset: float):
        self.form = form.lower()
        self.amplitude = amplitude
        self.frequency = frequency
        self.dc_offset = dc_offset
        self.omega = 2.0 * math.pi * self.frequency

    def __call__(self, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Расчет напряжения сигнала в момент(ы) времени t."""
        if self.form == 'sine':
            if isinstance(t, np.ndarray):
                return self.dc_offset + self.amplitude * np.sin(self.omega * t)
            return self.dc_offset + self.amplitude * math.sin(self.omega * float(t))
        elif self.form == 'square':
            # Меандр (ступенька): переключается между offset - amp и offset + amp
            # Обрабатываем как скаляры, так и массивы
            if isinstance(t, np.ndarray):
                sine_val = np.sin(self.omega * t)
                res = np.ones_like(t) * (self.dc_offset - self.amplitude)
                res[sine_val >= 0] = self.dc_offset + self.amplitude
                return res
            else:
                sine_val = math.sin(self.omega * float(t))
                return self.dc_offset + self.amplitude if sine_val >= 0 else self.dc_offset - self.amplitude
        elif self.form == 'triangle':
            # Треугольный сигнал
            if isinstance(t, np.ndarray):
                sine_val = np.sin(self.omega * t)
                return self.dc_offset + self.amplitude * (2.0 / np.pi) * np.arcsin(sine_val)
            sine_val = math.sin(self.omega * float(t))
            return self.dc_offset + self.amplitude * (2.0 / math.pi) * math.asin(sine_val)
        else:
            raise ValueError(f"Неизвестная форма сигнала: {self.form}")

CAPACITOR_PRESETS_PATH = Path(__file__).with_name("capacitors.json")

def load_capacitor_presets(path: Optional[Path] = None) -> Dict[str, Dict[str, object]]:
    """Load capacitor presets from the shared JSON catalog.

    The catalog uses SI units: capacitance is stored in farads.
    """
    catalog_path = path or CAPACITOR_PRESETS_PATH
    with catalog_path.open("r", encoding="utf-8") as f:
        return json.load(f)

CAPACITOR_PRESETS = load_capacitor_presets()

class CapacitorModel:
    """Модель нелинейной емкости конденсатора с интерполяцией LUT или аналитическим фитингом."""
    def __init__(self, cap_type: str, c_nom: float, v_rated: float = 10.0, custom_drop: Optional[float] = None, v50: Optional[float] = None, n: Optional[float] = None):
        self.cap_type = cap_type.upper()
        self.c_nom = c_nom
        self.v_rated = v_rated
        self.custom_drop = custom_drop
        self.v50 = v50
        self.n = n
        self._constant_c = self.c_nom if self.cap_type == 'C0G' else None

        if self.cap_type != 'PRESET_FIT':
            # Инициализация LUT таблиц (нормализованные по напряжению v_rated)
            # Для отрицательных напряжений характеристика симметрична
            v_lut_norm = np.array([-2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
            
            if self.cap_type == 'C0G':
                # Линейный конденсатор (емкость не зависит от напряжения)
                c_lut_norm = np.ones_like(v_lut_norm)
            elif self.cap_type == 'X7R':
                # Типичный профиль X7R: спад емкости до ~50% при номинальном напряжении
                c_lut_norm = np.array([0.15, 0.25, 0.50, 0.75, 0.90, 0.97, 1.0, 0.97, 0.90, 0.75, 0.50, 0.25, 0.15])
            elif self.cap_type == 'X5R':
                # Типичный профиль X5R: спад емкости до ~20% при номинальном напряжении
                c_lut_norm = np.array([0.05, 0.10, 0.20, 0.40, 0.60, 0.85, 1.0, 0.85, 0.60, 0.40, 0.20, 0.10, 0.05])
            elif self.cap_type == 'CUSTOM':
                # Пользовательский профиль: падение на custom_drop% при v_rated
                drop = custom_drop if custom_drop is not None else 0.50
                # Квадратичное падение емкости с ограничением снизу в 5% от номинала
                c_lut_norm = 1.0 - drop * (v_lut_norm ** 2)
                c_lut_norm = np.clip(c_lut_norm, 0.05, 1.0)
            else:
                raise ValueError(f"Неизвестный тип диэлектрика: {self.cap_type}")

            # Абсолютные значения напряжения для LUT
            self.v_points = v_lut_norm * self.v_rated
            self.c_points = c_lut_norm * self.c_nom
            self._v_points_list = self.v_points.tolist()
            self._c_points_list = self.c_points.tolist()

    def _interp_scalar(self, value: float) -> float:
        v_points = self._v_points_list
        c_points = self._c_points_list
        if value <= v_points[0]:
            return c_points[0]
        if value >= v_points[-1]:
            return c_points[-1]

        idx = bisect_left(v_points, value)
        v0 = v_points[idx - 1]
        v1 = v_points[idx]
        c0 = c_points[idx - 1]
        c1 = c_points[idx]
        return c0 + (c1 - c0) * ((value - v0) / (v1 - v0))

    def get_C(self, v: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Возвращает емкость конденсатора при напряжении v."""
        if self.cap_type == 'PRESET_FIT':
            v50 = self.v50 if self.v50 is not None and self.v50 > 0 else 5.0
            n = self.n if self.n is not None else 1.5
            if isinstance(v, np.ndarray):
                return self.c_nom / (1.0 + (np.abs(v) / v50) ** n)
            return self.c_nom / (1.0 + (abs(float(v)) / v50) ** n)
        else:
            if self._constant_c is not None:
                if isinstance(v, np.ndarray):
                    return np.full_like(v, self._constant_c, dtype=float)
                return self._constant_c
            if not isinstance(v, np.ndarray):
                return self._interp_scalar(float(v))
            return np.interp(v, self.v_points, self.c_points)

def _require_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")

def validate_simulation_parameters(
    R_ext: float,
    R_sw: float,
    C_sh: float,
    T_acq: float,
    f_s: float,
    N_samples: int,
    V_ref: float,
    steps_per_cycle: int
) -> None:
    """Validate physical and numerical limits before running a simulation."""
    _require_positive("R_ext", R_ext)
    _require_positive("R_sw", R_sw)
    _require_positive("C_sh", C_sh)
    _require_positive("T_acq", T_acq)
    _require_positive("f_s", f_s)
    _require_positive("V_ref", V_ref)

    if N_samples <= 0:
        raise ValueError(f"N_samples must be positive, got {N_samples}")
    if steps_per_cycle <= 1:
        raise ValueError(f"steps_per_cycle must be greater than 1, got {steps_per_cycle}")

    T_s = 1.0 / f_s
    if T_acq >= T_s:
        raise ValueError(f"T_acq must be shorter than the sampling period ({T_s:g} s), got {T_acq:g} s")

def _advance_closed_linearized(
    v_ext: float,
    v_sh: float,
    v_in: float,
    dt: float,
    r_ext: float,
    r_sw: float,
    c_ext: float,
    c_sh: float
) -> Tuple[float, float]:
    """Exact 2x2 RC step with C_ext frozen at the start of the step."""
    a = 1.0 / (r_ext * c_ext)
    b = 1.0 / (r_sw * c_ext)
    c = 1.0 / (r_sw * c_sh)
    trace = -(a + b + c)
    discriminant = max(0.0, trace * trace - 4.0 * a * c)
    eigen_gap = math.sqrt(discriminant)
    lambda_1 = 0.5 * (trace + eigen_gap)
    lambda_2 = 0.5 * (trace - eigen_gap)
    exp_1 = math.exp(lambda_1 * dt) if lambda_1 * dt > -745.0 else 0.0
    exp_2 = math.exp(lambda_2 * dt) if lambda_2 * dt > -745.0 else 0.0

    delta_ext = v_ext - v_in
    delta_sh = v_sh - v_in
    if eigen_gap <= 1e-18 * max(1.0, abs(trace)):
        eigenvalue = 0.5 * trace
        decay = math.exp(eigenvalue * dt) if eigenvalue * dt > -745.0 else 0.0
        next_ext = decay * (
            delta_ext + dt * ((-a - b - eigenvalue) * delta_ext + b * delta_sh)
        )
        next_sh = decay * (
            delta_sh + dt * (c * delta_ext + (-c - eigenvalue) * delta_sh)
        )
    else:
        alpha = (exp_1 - exp_2) / eigen_gap
        beta = (-lambda_2 * exp_1 + lambda_1 * exp_2) / eigen_gap
        next_ext = alpha * ((-a - b) * delta_ext + b * delta_sh) + beta * delta_ext
        next_sh = alpha * (c * delta_ext - c * delta_sh) + beta * delta_sh

    return v_in + next_ext, v_in + next_sh

def simulate_adc_input(
    signal_gen: SignalGenerator,
    cap_model: CapacitorModel,
    R_ext: float,
    R_sw: float,
    C_sh: float,
    T_acq: float,
    f_s: float,
    N_samples: int,
    V_ref: float,
    V_ext_init: float = 0.0,
    V_sh_init: float = 0.0,
    method: str = 'semi_implicit',
    steps_per_cycle: int = 100,
    store_trace: bool = True,
    trace_until: Optional[float] = None
) -> Dict[str, np.ndarray]:
    """
    Симуляция входного тракта АЦП во временной области.
    
    Параметры:
        signal_gen: объект SignalGenerator
        cap_model: объект CapacitorModel для C_ext
        R_ext: внешнее сопротивление (Ом)
        R_sw: сопротивление ключа выборки (Ом)
        C_sh: емкость выборки АЦП (Ф)
        T_acq: время выборки (с)
        f_s: частота дискретизации АЦП (Гц)
        N_samples: число выборок
        V_ref: опорное напряжение АЦП (В)
        V_ext_init: начальное напряжение на C_ext (В)
        V_sh_init: начальное напряжение на C_sh (В)
        method: метод решения ('semi_implicit' или 'radau')
        steps_per_cycle: число шагов на один цикл выборки (для semi_implicit)
        store_trace: сохранять ли полные временные массивы для графиков
        trace_until: опциональная верхняя граница времени для сохранения траектории
        
    Возвращает:
        Словарь с массивами результатов:
            't': время непрерывного процесса
            'v_in': входной сигнал
            'v_ext': напряжение на внешнем конденсаторе C_ext
            'v_sh': напряжение на внутреннем конденсаторе C_sh
            'is_closed': состояние ключа (1 - замкнут, 0 - разомкнут)
            't_samples': моменты фиксации выборок АЦП
            'v_samples': аналоговые напряжения в моменты фиксации
            'digital_codes': оцифрованные N-битные коды (клиппированные)
    """
    validate_simulation_parameters(R_ext, R_sw, C_sh, T_acq, f_s, N_samples, V_ref, steps_per_cycle)
    T_s = 1.0 / f_s
    
    if method == 'semi_implicit':
        # Для полунеявной схемы Эйлера разобьем цикл на фиксированные шаги
        # Выделим целое число шагов для фазы выборки и фазы преобразования
        N_acq = max(3, int(np.round(steps_per_cycle * (T_acq / T_s))))
        N_conv = max(1, steps_per_cycle - N_acq)
        
        dt_acq = T_acq / N_acq
        dt_conv = (T_s - T_acq) / N_conv
        
        # Резервируем массивы под траектории только когда они нужны для графиков.
        if store_trace:
            trace_cycles = N_samples
            if trace_until is not None:
                trace_cycles = min(N_samples, max(0, int(math.ceil(trace_until / T_s))))
            total_steps = trace_cycles * (N_acq + N_conv)
            t_trace = np.zeros(total_steps)
            v_in_trace = np.zeros(total_steps)
            v_ext_trace = np.zeros(total_steps)
            v_sh_trace = np.zeros(total_steps)
            is_closed_trace = np.zeros(total_steps)
        else:
            t_trace = np.empty(0)
            v_in_trace = np.empty(0)
            v_ext_trace = np.empty(0)
            v_sh_trace = np.empty(0)
            is_closed_trace = np.empty(0)
        
        t_samples = np.zeros(N_samples)
        v_samples = np.zeros(N_samples)
        v_ext_samples = np.zeros(N_samples)
        
        V_ext, V_sh = V_ext_init, V_sh_init
        G_ext = 1.0 / R_ext
        
        idx = 0
        for n in range(N_samples):
            t_cycle_start = n * T_s
            
            # --- 1. Фаза выборки (Switch Closed) ---
            a_sh = dt_acq / (R_sw * C_sh)
            G_sw = 1.0 / (R_sw * (1.0 + a_sh))
            for step in range(N_acq):
                t_curr = t_cycle_start + (step + 1) * dt_acq
                t_step_start = t_curr - dt_acq
                v_in_val = signal_gen(t_curr)
                
                # Полунеявный шаг: используем C_ext(V_ext) с предыдущего шага
                C_e = max(cap_model.get_C(V_ext), 1e-15)

                # Moderate external time constants benefit from the exact linearized
                # coupled step. Very stiff nonlinear cases retain L-stable backward Euler.
                if dt_acq / (R_ext * C_e) < 10.0:
                    v_in_mid = signal_gen(t_step_start + 0.5 * dt_acq)
                    V_ext, V_sh = _advance_closed_linearized(
                        V_ext, V_sh, v_in_mid, dt_acq, R_ext, R_sw, C_e, C_sh
                    )
                else:
                    num = V_ext + (dt_acq / C_e) * (G_ext * v_in_val + G_sw * V_sh)
                    den = 1.0 + (dt_acq / C_e) * (G_ext + G_sw)
                    V_ext = num / den
                    V_sh = (V_sh + a_sh * V_ext) / (1.0 + a_sh)
                
                if store_trace and idx < len(t_trace):
                    t_trace[idx] = t_curr
                    v_in_trace[idx] = v_in_val
                    v_ext_trace[idx] = V_ext
                    v_sh_trace[idx] = V_sh
                    is_closed_trace[idx] = 1.0
                    idx += 1
            
            # Точка фиксации выборки АЦП (конец фазы выборки)
            t_samples[n] = t_cycle_start + T_acq
            v_samples[n] = V_sh
            v_ext_samples[n] = V_ext
            
            # --- 2. Фаза преобразования (Switch Open) ---
            for step in range(N_conv):
                t_curr = t_cycle_start + T_acq + (step + 1) * dt_conv
                v_in_val = signal_gen(t_curr)
                
                C_e = max(cap_model.get_C(V_ext), 1e-15)
                
                num = V_ext + (dt_conv / C_e) * G_ext * v_in_val
                den = 1.0 + (dt_conv / C_e) * G_ext
                V_ext = num / den
                # V_sh не меняется (разомкнуто)
                
                if store_trace and idx < len(t_trace):
                    t_trace[idx] = t_curr
                    v_in_trace[idx] = v_in_val
                    v_ext_trace[idx] = V_ext
                    v_sh_trace[idx] = V_sh
                    is_closed_trace[idx] = 0.0
                    idx += 1
                
        return {
            't': t_trace,
            'v_in': v_in_trace,
            'v_ext': v_ext_trace,
            'v_sh': v_sh_trace,
            'is_closed': is_closed_trace,
            't_samples': t_samples,
            'v_samples': v_samples,
            'v_ext_samples': v_ext_samples
        }
        
    elif method == 'radau':
        # Поцикловое решение через жесткий ОДУ решатель Radau (scipy)
        t_trace_list = []
        v_in_trace_list = []
        v_ext_trace_list = []
        v_sh_trace_list = []
        is_closed_list = []
        
        t_samples = np.zeros(N_samples)
        v_samples = np.zeros(N_samples)
        v_ext_samples = np.zeros(N_samples)
        
        V_ext, V_sh = V_ext_init, V_sh_init
        
        for n in range(N_samples):
            t_cycle_start = n * T_s
            t_acq_time = t_cycle_start + T_acq
            t_cycle_end = (n + 1) * T_s
            
            # --- 1. Фаза выборки (Switch Closed) ---
            def ode_closed(t, y):
                V_e, V_s = y
                C_e = max(cap_model.get_C(V_e), 1e-15)
                v_in_val = signal_gen(t)
                dV_e = ( (v_in_val - V_e)/R_ext - (V_e - V_s)/R_sw ) / C_e
                dV_s = (V_e - V_s) / (R_sw * C_sh)
                return [dV_e, dV_s]
            
            # Генерируем равномерные точки для детального графика внутри фазы
            t_eval_acq = np.linspace(t_cycle_start, t_acq_time, 15)
            sol1 = integrate.solve_ivp(
                ode_closed, 
                [t_cycle_start, t_acq_time], 
                [V_ext, V_sh], 
                method='Radau', 
                t_eval=t_eval_acq,
                rtol=1e-6, 
                atol=1e-8
            )
            
            V_ext, V_sh = sol1.y[0][-1], sol1.y[1][-1]
            
            # Запоминаем траектории
            t_trace_list.append(sol1.t)
            v_ext_trace_list.append(sol1.y[0])
            v_sh_trace_list.append(sol1.y[1])
            v_in_trace_list.append(signal_gen(sol1.t))
            is_closed_list.append(np.ones_like(sol1.t))
            
            t_samples[n] = t_acq_time
            v_samples[n] = V_sh
            v_ext_samples[n] = V_ext
            
            # --- 2. Фаза преобразования (Switch Open) ---
            def ode_open(t, y):
                V_e, V_s = y
                C_e = max(cap_model.get_C(V_e), 1e-15)
                v_in_val = signal_gen(t)
                dV_e = (v_in_val - V_e) / (R_ext * C_e)
                dV_s = 0.0
                return [dV_e, dV_s]
                
            t_eval_conv = np.linspace(t_acq_time, t_cycle_end, 35)
            sol2 = integrate.solve_ivp(
                ode_open, 
                [t_acq_time, t_cycle_end], 
                [V_ext, V_sh], 
                method='Radau', 
                t_eval=t_eval_conv,
                rtol=1e-6, 
                atol=1e-8
            )
            
            V_ext, V_sh = sol2.y[0][-1], sol2.y[1][-1]
            
            t_trace_list.append(sol2.t)
            v_ext_trace_list.append(sol2.y[0])
            v_sh_trace_list.append(sol2.y[1])
            v_in_trace_list.append(signal_gen(sol2.t))
            is_closed_list.append(np.zeros_like(sol2.t))
            
        return {
            't': np.concatenate(t_trace_list),
            'v_in': np.concatenate(v_in_trace_list),
            'v_ext': np.concatenate(v_ext_trace_list),
            'v_sh': np.concatenate(v_sh_trace_list),
            'is_closed': np.concatenate(is_closed_list),
            't_samples': t_samples,
            'v_samples': v_samples,
            'v_ext_samples': v_ext_samples
        }
    else:
        raise ValueError(f"Неизвестный метод симуляции: {method}")

def quantize_samples(v_samples: np.ndarray, V_ref: float, bits: int) -> Tuple[np.ndarray, bool]:
    """
    Квантование и ограничение (clipping) аналогового сигнала в N-битный код.
    
    Параметры:
        v_samples: массив аналоговых отсчетов (В)
        V_ref: опорное напряжение АЦП (В)
        bits: разрядность АЦП (целое число, например, 12, 16, 24)
        
    Возвращает:
        (codes, clipping_detected):
            codes: массив целых чисел от 0 до 2^bits - 1
            clipping_detected: True, если обнаружено ограничение (выход за пределы 0 или V_ref)
    """
    _require_positive("V_ref", V_ref)
    if bits <= 0:
        raise ValueError(f"bits must be positive, got {bits}")

    # Проверка на клиппинг (выход за пределы 0 и V_ref)
    clipping_detected = np.any(v_samples >= V_ref) or np.any(v_samples <= 0.0)
    
    # Жесткое ограничение
    v_clipped = np.clip(v_samples, 0.0, V_ref)
    
    # Квантование
    levels = (2 ** bits) - 1
    codes = np.round((v_clipped / V_ref) * levels).astype(np.int64)
    
    return codes, clipping_detected

def reconstruct_voltage(codes: np.ndarray, V_ref: float, bits: int) -> np.ndarray:
    """Восстановление напряжения из цифровых кодов АЦП."""
    _require_positive("V_ref", V_ref)
    if bits <= 0:
        raise ValueError(f"bits must be positive, got {bits}")

    levels = (2 ** bits) - 1
    return codes.astype(np.float64) * (V_ref / levels)
