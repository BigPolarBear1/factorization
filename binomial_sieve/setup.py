import sys

from setuptools import Extension,setup
from Cython.Build import cythonize

ext = Extension("binsieve",["binsieve.pyx"],include_dirs=sys.path,extra_compile_args=['-O3','-march=native'])  #libraries=['gmp','mpfr','mpc']

ext.cython_directives={'language_level':"3",'profile':True,'linetrace':True}

setup(name="binsieve",ext_modules=cythonize([ext],include_path=sys.path))