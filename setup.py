from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "monod_cpp",
        ["monod_bindings.cpp", "monod_model.cpp"],
        cxx_std=11,
    ),
]

setup(
    name="monod-pybind11",
    version="0.1.0",
    description="Pybind11 bindings for Monod bioreaction model",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
