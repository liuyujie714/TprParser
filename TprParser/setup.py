from setuptools import setup, Extension
import numpy, os

# Always update program version
__version__ = '0.11.7'


# set c++17
CXXFLAGS = []
if os.name == 'nt':
    CXXFLAGS = ['/std:c++17', '/O2']
elif os.name == 'posix':
    os.environ['CC'] = 'g++' # tell compiler cpp must need g++
    CXXFLAGS = ['-std=c++17', '-O3']
else:
    raise OSError("Unsupported OS %s" % os.name)

module = Extension(
    name='TprParser',  # module name
    include_dirs=[numpy.get_include()], # need numpy
    language='C++',
    extra_compile_args=CXXFLAGS, # C++ standard
    define_macros=[('_CRT_SECURE_NO_WARNINGS', 1)], # for MSVC
    sources=['Py_Module.cpp', 'Reader.cpp'], # source code path
)
setup(
    name='TprParser',
    version=__version__,
    description='A reader of gromacs tpr file', 
    author='Yujie Liu',
    author_email='',
    python_requires='>=3.8',
    install_requires=['typing_extensions', 'numpy'],
    packages=['.'],
    exclude=['setup.py'],
    ext_modules=[module]
)
