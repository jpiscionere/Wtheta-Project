from distutils.core import setup, Extension
import numpy.distutils.misc_util


from distutils.command.build_ext import build_ext
copt =  {'msvc': ['/openmp', '/Ox', '/fp:fast','/favor:INTEL64','/Og']  ,
         'mingw32' : ['-fopenmp','-O3','-ffast-math','-march=native'],
         'unix': ['-std=c99','-Wall','-Wextra','-Wshadow']}
lopt =  {'mingw32' : ['-fopenmp'],
         'unix': ['']}

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        c = self.compiler.compiler_type
        if copt.has_key(c):
            for e in self.extensions:
                e.extra_compile_args = copt[ c ]
        if lopt.has_key(c):
             for e in self.extensions:
                 e.extra_link_args = lopt[ c ]
        build_ext.build_extensions(self)



setup(
    ext_modules=[Extension("_chi2_fof", ["_chi2_fof.c", "chi2_fof.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)


