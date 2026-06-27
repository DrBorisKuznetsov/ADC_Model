# ADC Input Stage Simulator (ADC Input Model)

[![Python Tests](https://github.com/DrBorisKuznetsov/ADC_Model/actions/workflows/python-tests.yml/badge.svg?branch=main)](https://github.com/DrBorisKuznetsov/ADC_Model/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An interactive simulation tool for the Analog-to-Digital Converter (ADC) input stage. It models the physical effects of non-linear filter capacitors (such as the DC bias effect in MLCC Class II dielectrics like X7R and X5R) and sampling switch settling dynamics on ADC performance metrics (THD, SNR, SINAD, ENOB).

The simulator is built with a Python mathematical core running in the browser via WebAssembly (PyScript/Pyodide) for zero-dependency client-side execution.

---

## Features

1. **Time-Domain Simulation**:
   - Accurately models the sample-and-hold (S&H) circuit behavior during the *acquisition* (switch closed) and *conversion* (switch open) phases.
   - Utilizes stiff ODE solvers: high-accuracy implicit `Radau` solver (`scipy`) and a custom, ultra-fast semi-implicit scheme with exact 2x2 local linearization.
2. **Non-linear Capacitor Model ($C(V)$)**:
   - Supports linear `C0G` and non-linear `X7R`, `X5R`, `CUSTOM`, and `PRESET_FIT` capacitor profiles.
   - Includes real commercial MLCC presets (Murata, TDK, Samsung, AVX, Kemet, Taiyo Yuden) loaded dynamically from a unified catalog.
3. **Spectral Analysis (FFT)**:
   - Computes THD, SNR, SINAD, and ENOB using a 4-term Blackman-Harris window with noise power bandwidth correction.
4. **ENOB Loss Diagnostics & Ablation**:
   - Decomposes total dynamic range loss into specific contributors: $C(V)$ non-linearity, acquisition settling error, and their cross-coupling interaction.
5. **N-Dimensional Design Space Sweep**:
   - Generates interactive heatmaps to map safe operating regions across parameter combinations (e.g., $R_{ext}$ vs. $C_{ext}$).
6. **PDF Reports**:
   - Generates and downloads vector PDF reports including all numeric results, active settings, and Plotly chart snapshots.

---

## Web Application

The interactive simulator is designed to run as a client-side Single-Page Application (SPA) and is published online as a web resource. Since all computations (including the Python-based solver) execute inside the browser using WebAssembly, the client browser handles 100% of the workload.

For local offline development and report testing, the workspace includes a helper server script:
```bash
python serve_app.py
```

---

## Python Math Core

The underlying simulation algorithms can also be run locally using the Python source files for batch simulations and verification.

### Installation
Install the required dependencies:
```bash
pip install -r requirements.txt
```

### Running Unit Tests
Validate the solver, quantizer, and spectral analysis components:
```bash
python -m unittest discover -s tests
```

### Physical Verification
Run the verification script to simulate and compare C0G, X7R, and X5R dielectrics under identical conditions:
```bash
python -m src.verify
```
The results and charts are saved to `wiki/assets/verification_results.png`.

---

## Documentation

Comprehensive mathematical formulations, physical schemas, and component details are maintained in the project's LLM-Wiki:
* **Index**: [wiki/index.md](file:///g:/%D0%9C%D0%BE%D0%B9%20%D0%B4%D0%B8%D1%81%D0%BA/SNR_Lib/ADC_Model/wiki/index.md)
* **Log**: [wiki/log.md](file:///g:/%D0%9C%D0%BE%D0%B9%20%D0%B4%D0%B8%D1%81%D0%BA/SNR_Lib/ADC_Model/wiki/log.md)

---

## Channel Reference
Created in collaboration with the engineering YouTube channel [@High_SNR_Channel](https://www.youtube.com/@High_SNR_Channel).
