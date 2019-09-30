try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.command.sdist import sdist as _sdist
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from distutils.extension import Extension

import numpy as np

setup(
    name='pyflagstats',
    version='0.1.0',
    description="Efficient subroutines for computing summary statistics for the SAM FLAG field",
    long_description="""
Given a stream of k$bit words, we seek to sum the bit values at indexes 0, 1, 2, ..., k-1 across multiple words 
by computing k distinct sums. If the k-bit words are one-hot encoded then the sums corresponds to their frequencies.

This multiple-sum problem is a generalization of the population-count 
problem where we count the total number of set bits in independent machine words.
We refer to this new problem as the positional population-count problem.

Using SIMD (Single Instruction, Multiple Data) instructions from 
recent Intel processors, we describe algorithms for computing the 16-bit position population count using about one eighth (0.125) 
of a CPU cycle per 16-bit word. Our best approach is about 140-fold faster than competitive code using only non-SIMD instructions in terms of CPU cycles.

This package contains the application of the efficient positional population count operator to computing summary statistics for the SAM FLAG field.""",
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
    keywords = ['simd', 'popcount', 'popcnt', 'pospopcnt', 'hts', 'ngs', 'flags']
)