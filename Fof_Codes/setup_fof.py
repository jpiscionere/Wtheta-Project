from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    ext_modules=[Extension("_chi2_fof", ["_chi2_fof.c", "chi2_fof.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)


