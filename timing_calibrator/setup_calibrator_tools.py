#!/usr/bin/env python3
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext = Extension("calibrator_tools", ["calibrator_tools.pyx"],
    include_dirs=[np.get_include()],
    library_dirs=[],
    libraries=[]
)
 
setup(ext_modules=[ext],
    cmdclass = {'build_ext': build_ext})

# from LoLIM.utilities import GSL_include, GSL_library_dir

# ext = Extension("calibrator_tools", ["calibrator_tools.pyx"],
#     include_dirs=[np.get_include(), 
#                   GSL_include()],
#     library_dirs=[GSL_library_dir()],
#     libraries=["gsl", 'blas']
# )
 
# setup(ext_modules=[ext],
#     cmdclass = {'build_ext': build_ext})