# Monod bioreaction (Python + C++11 + pybind11)

This repository contains an implementation of the Monod bioreaction model in Python and C++11, with a pybind11 interface to call the C++ core from Python.

Repository URL: `https://github.com/michelanthony/TestAICoding.git`

![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]
(https://colab.research.google.com/github/michelanthony/TestAICoding/notebook/notebook.ipynb)

## Folder structure

- `monod_bioreaction.py`
  Main Python script: Monod equations, Euler integration in Python, fallback when the C++ module is not available, and tabulated output generation.

- `monod_model.hpp` / `monod_model.cpp`
  C++11 model core (`dX_dt`, `dS_dt`, `integrate`, and `SimulationResult`).

- `monod_bindings.cpp`
  pybind11 bindings exposing the C++ core as Python module `monod_cpp`.

- `setup.py` / `pyproject.toml`
  Build files to compile the pybind11 extension (`monod_cpp`) with C++11.

- `test_monod_cpp_main.cpp`
  Standalone C++ test (`main`) that compares C++ results with Python reference values and writes a tabulated file.

- `generate_python_reference.py`
  Utility script to regenerate Python reference values used by the C++ test.

- `notebooks/monod_colab_tutorial.ipynb`
  Google Colab tutorial notebook that generates plots with Matplotlib and saves them into `results/`.

## Run and test

### 1) C++ test (numeric comparison)

```bash
g++ -std=c++11 -Wall -Wextra -pedantic monod_model.cpp test_monod_cpp_main.cpp -o test_monod_cpp_main
./test_monod_cpp_main
```

The executable:

- prints a detailed comparison table in the terminal,
- creates a `results` folder,
- writes `results/monod_test_results.tsv` with numeric values,
- and returns success if vectors `time`, `biomass`, and `substrate` match the Python reference.

### 2) Build pybind11 extension

```bash
pip install -e .
```

Requires `pybind11`, `setuptools`, and a C++11 compiler.

### 3) Run Python simulation (tabulated output)

```bash
python monod_bioreaction.py
```

The script uses `monod_cpp` when available, otherwise it uses the pure Python implementation.
It writes `results/monod_python_results.tsv`.
The base script does not require NumPy for execution.

## Google Colab tutorial

A Colab-ready tutorial notebook is available at:

- `notebooks/monod_colab_tutorial.ipynb`

What it includes:

- install and setup steps for Colab,
- Python and pybind11 simulation calls,
- tabulated output generation in `results/monod_notebook_comparison.tsv`,
- Matplotlib plot generation in `results/monod_notebook_plot.png`.

The notebook is the place where plots are generated and added to the `results` folder.

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
git clone https://github.com/michelanthony/TestAICoding.git
cd TestAICoding
git add .
git commit -m "your change summary"
git push origin <your-branch>
```

Then open a pull request on GitHub from `<your-branch>` to your target branch.
