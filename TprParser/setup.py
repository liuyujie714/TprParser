from setuptools import setup, Extension
import numpy, os

# Always update program version
__version__ = '0.1.41'

# Description
long_doc = ""
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_doc = f.read()
except:
    pass
    

# set c++17
CXXFLAGS = []
LINKER   = []
if os.name == 'nt':
    CXXFLAGS = ['/std:c++17', '/O2']
elif os.name == 'posix':
    os.environ['CC'] = 'g++' # tell compiler cpp must need g++
    CXXFLAGS = ['-std=c++17', '-O3']
    LINKER = ['-lstdc++'] 
else:
    raise OSError("Unsupported OS %s" % os.name)

module = Extension(
    name='TprParser_',  # module name from PyMODINIT_FUNC PyInit_TprParser_(void) 
    include_dirs=[numpy.get_include()], # need numpy
    language='c++',
    extra_compile_args=CXXFLAGS, # C++ standard
    define_macros=[('_CRT_SECURE_NO_WARNINGS', 1)], # for MSVC
    sources=['src/Py_Module.cpp', 'src/Reader.cpp'], # source code path
    extra_link_args=LINKER  # link to c++ library
)
setup(
    name='TprParser',
    version=__version__,
    description='A reader of gromacs tpr file', 
    long_description=long_doc,
    long_description_content_type='text/markdown',
    license='',
    author='Yujie Liu',
    author_email='',
    python_requires='>=3.8',
    install_requires=['typing_extensions', 'numpy'],
    exclude=['setup.py'],
    ext_modules=[module],
    # put TprReader.py/__init__.py in TprParser folder to site-packages
    py_modules=['TprParser.TprReader', 'TprParser.__init__'] 
)
