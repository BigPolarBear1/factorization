import sys

from setuptools import Extension,setup
from Cython.Build import cythonize

ext = Extension("QSv3_simd",["QSv3_simd.pyx"],include_dirs=sys.path,extra_compile_args=['-O3','-march=native'])  #libraries=['gmp','mpfr','mpc']

ext.cython_directives={'language_level':"3",'profile':True,'linetrace':True}

setup(name="QSv3_simd",ext_modules=cythonize([ext],include_path=sys.path))