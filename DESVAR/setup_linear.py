from distutils.core import setup, Extension
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("lnlike_fast_linear", ["lnlike_fast_linear.pyx"])]
)
