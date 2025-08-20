from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("metrics", ["metrics.pyx"],
    include_dirs=[numpy.get_include()])]
)
