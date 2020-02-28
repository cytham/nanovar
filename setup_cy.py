from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize('nanovar/nv_bam_parser.pyx'))
