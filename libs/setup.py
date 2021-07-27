from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
    ext_modules = cythonize("distances.pyx"),
    include_dirs=[numpy.get_include()]
)

#Run: python3 setup.py build_ext --inplace
