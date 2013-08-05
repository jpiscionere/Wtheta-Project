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
                        
                        
module=Extension("_chi2",
                  include_dirs=['../utils'],
                  define_macros = [('MAJOR_VERSION', '1'),
                                   ('MINOR_VERSION', '0')],
                  sources=["_chi2.c", "chi2.c","../utils/utils.c","../utils/ftread.c","../utils/bgc_read_utils.c"]
                  
                  )

setup(
    name="LasDamas wp Fitting routines",
    version = '1.0',
    description = 'Fits SDSS data to LasDamas boxes',
    author = 'Manodeep Sinha',
    author_email = 'manodeep@gmail.com',
    url = 'http://astro.phy.vanderbilt.edu/~sinham/',
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    ext_modules=[module],
    cmdclass = {'build_ext': build_ext_subclass },
)


