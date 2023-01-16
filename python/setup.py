from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


icsi2023_extension = Extension(
    name="icsi2023",
    sources=["icsi2023.pyx"],
    libraries=["icsi2023"],
    library_dirs=["lib"],
    include_dirs=["lib"],
)
setup(name="icsi2023", ext_modules=cythonize([icsi2023_extension]), include_dirs=[numpy.get_include()], compiler_directives={'language_level' : "3"})
