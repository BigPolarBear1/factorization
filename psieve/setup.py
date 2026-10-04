import sys

from setuptools import Extension,setup
from Cython.Build import cythonize

args=['-O3','-march=native']
exts = [Extension("QSv3_simd",["QSv3_simd.pyx"],include_dirs=sys.path,extra_compile_args=args),  #libraries=['gmp','mpfr','mpc']
        Extension("bitsieve",["bitsieve.pyx"],extra_compile_args=args)]

for ext in exts:
    ext.cython_directives={'language_level':"3",'profile':False,'linetrace':False}

setup(name="QSv3_simd",ext_modules=cythonize(exts,include_path=sys.path))