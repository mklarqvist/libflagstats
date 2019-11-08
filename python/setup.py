try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.command.sdist import sdist as _sdist
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from distutils.extension import Extension

import numpy as np

# Read contents of README markdown file and store
# in the long_description parameter
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='UTF-8') as f:
    long_description = f.read()

setup(
    name='pyflagstats',
    version='0.1.4',
    description="Efficient subroutines for computing summary statistics for the SAM FLAG field",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Marcus D. R. Klarqvist",
    author_email="mk819@cam.ac.uk",
    platforms="Linux, MacOSX, Windows",
    url="https://github.com/mklarqvist/libflagstats",
    ext_modules=cythonize(
        Extension(
            "pyflagstats",
            sources=["libflagstats.pyx"],
            include_dirs=[np.get_include()]
        )
    ),
    install_requires=["numpy", "cython"],
    license="Apache 2.0",
    keywords = ['simd', 'popcount', 'popcnt', 'pospopcnt', 'hts', 'ngs', 'flags'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
)
