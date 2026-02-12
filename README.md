# Monod bioreaction (Python + C++11 + pybind11)

This repository contains an implementation of the Monod bioreaction model in Python and C++11, with a pybind11 interface to call the C++ core from Python.

## Folder structure

- `monod_bioreaction.py`
  Main Python script: Monod equations, Euler integration in Python, fallback when the C++ module is not available, and curve plotting.

- `monod_model.hpp` / `monod_model.cpp`
  C++11 model core (`dX_dt`, `dS_dt`, `integrate`, and `SimulationResult`).

- `monod_bindings.cpp`
  pybind11 bindings exposing the C++ core as Python module `monod_cpp`.

- `setup.py` / `pyproject.toml`
  Build files to compile the pybind11 extension (`monod_cpp`) with C++11.

- `test_monod_cpp_main.cpp`
  Standalone C++ test (`main`) that compares C++ results with Python reference values.

- `generate_python_reference.py`
  Utility script to regenerate Python reference values used by the C++ test.

## Run and test

### 1) C++ test (numeric comparison)

```bash
g++ -std=c++11 -Wall -Wextra -pedantic monod_model.cpp test_monod_cpp_main.cpp -o test_monod_cpp_main
./test_monod_cpp_main
```

The executable:

- prints a detailed comparison table in the terminal,
- creates a `results` folder,
- writes `results/monod_test_results.csv` with numeric values,
- writes `results/monod_test_plot.svg` with overlay plots (Python reference vs C++),
- and returns success if vectors `time`, `biomass`, and `substrate` match the Python reference.

### 2) Build pybind11 extension

```bash
pip install -e .
```

Requires `pybind11`, `setuptools`, and a C++11 compiler.

### 3) Run Python simulation

```bash
python monod_bioreaction.py
```

The script uses `monod_cpp` when available, otherwise it uses the pure Python implementation.

## Curve comparison (Python vs C++)

The snippet below overlays curves from Python integration and C++ integration:

```python
import numpy as np
import matplotlib.pyplot as plt
import monod_bioreaction as m

X0, S0, t_max, dt = 0.1, 1.0, 10.0, 0.1

t_py, x_py, s_py = m.integrate(X0, S0, t_max, dt)

if m.monod_cpp is None:
    raise RuntimeError("monod_cpp module is not available. Run `pip install -e .`.")

t_cpp, x_cpp, s_cpp = m.integrate_with_pybind11(X0, S0, t_max, dt)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(t_py, x_py, label="Biomass Python", linestyle="--")
plt.plot(t_cpp, x_cpp, label="Biomass C++", alpha=0.8)
plt.xlabel("Time (h)")
plt.ylabel("Concentration (g/L)")
plt.title("Biomass comparison")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t_py, s_py, label="Substrate Python", linestyle="--")
plt.plot(t_cpp, s_cpp, label="Substrate C++", alpha=0.8)
plt.xlabel("Time (h)")
plt.ylabel("Concentration (g/L)")
plt.title("Substrate comparison")
plt.legend()

plt.tight_layout()
plt.show()
```

If implementations are consistent, Python and C++ curves overlap.


## CI/CD build artifacts (Linux and Windows)

GitHub Actions workflow `.github/workflows/build-artifacts.yml` builds the C++11 target on Linux and Windows.

It generates two artifacts:

- `monod-build-linux`
  - `test_monod_cpp_main`
  - `libmonod_model.a`

- `monod-build-windows`
  - `test_monod_cpp_main.exe`
  - `monod_model.lib`

The workflow also runs the C++ numeric comparison test before uploading artifacts.

## GitHub pull request

Use these commands after local validation:

```bash
git add .
git commit -m "your change summary"
git push origin <your-branch>
```

Then open a pull request on GitHub from `<your-branch>` to your target branch.
