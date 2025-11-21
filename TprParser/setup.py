from setuptools import setup, Extension
import numpy, os, re


# Always update program version
__version__ = '0.1.59'

# fix version
finit = 'TprParser/version.py'
with open(finit, 'r+') as f:
    lines = f.read()
    context = re.sub(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]', 
                     r'__version__ = "{}"'.format(__version__), lines)
    f.seek(0)
    f.write(context)
    f.truncate()


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

tpr_module = Extension(
    name='TprParser_',  # module name from PyMODINIT_FUNC PyInit_TprParser_(void) 
    include_dirs=[numpy.get_include()], # need numpy
    language='c++',
    extra_compile_args=CXXFLAGS, # C++ standard
    define_macros=[('_CRT_SECURE_NO_WARNINGS', 1)], # for MSVC
    sources=['src/Py_Tpr.cpp', 'src/Reader.cpp', 'src/Utils.cpp'], # source code path
    extra_link_args=LINKER  # link to c++ library
)

edr_module = Extension(
    name='EdrParser_',  # module name from PyMODINIT_FUNC PyInit_EdrParser_(void) 
    include_dirs=[numpy.get_include()], # need numpy
    language='c++',
    extra_compile_args=CXXFLAGS, # C++ standard
    define_macros=[('_CRT_SECURE_NO_WARNINGS', 1)], # for MSVC
    sources=['src/Py_Edr.cpp', 'src/EdrReader.cpp'], # source code path
    extra_link_args=LINKER  # link to c++ library
)

setup(
    name='TprParser',
    version=__version__,
    description='A reader of gromacs tpr file', 
    long_description=long_doc,
    long_description_content_type='text/markdown',
    license='GPL',
    author='Yujie Liu',
    author_email='',
    python_requires='>=3.8',
    install_requires=['typing_extensions<=4.12.2', 'numpy'],
    ext_modules=[tpr_module, edr_module],
    # put TprReader.py/__init__.py in TprParser folder to site-packages
    py_modules=['TprParser.TprReader', 
                'TprParser.__init__', 
                'TprParser.TprMakeTop',
                'TprParser.EdrReader',
                'TprParser.version'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
    ]
)
