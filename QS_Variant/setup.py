import sys

from setuptools import Extension,setup
from Cython.Build import cythonize

ext = Extension("QSv3_simd",["QSv3_simd.pyx"],include_dirs=sys.path,libraries=['gmp','mpfr','mpc'],extra_compile_args=['-O3','-march=native'])

ext.cython_directives={'language_level':"3"}

setup(name="QSv3_simd",ext_modules=cythonize([ext],include_path=sys.path))
