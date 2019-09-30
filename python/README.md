# pyflagstats

[![PyPI version](https://badge.fury.io/py/pyflagstats.svg)](https://badge.fury.io/py/pyflagstats)

Given a stream of k-bit words, we seek to sum the bit values at indexes 0, 1, 2,
..., k-1 across multiple words by computing k distinct sums. If the k-bit words
are one-hot encoded then the sums corresponds to their frequencies.

This multiple-sum problem is a generalization of the population-count problem
where we count the total number of set bits in independent machine words. We
refer to this new problem as the positional population-count problem.

Using SIMD (Single Instruction, Multiple Data) instructions from recent Intel
processors, we describe algorithms for computing the 16-bit position population
count using about one eighth (0.125) of a CPU cycle per 16-bit word. Our best
approach is about 140-fold faster than competitive code using only non-SIMD
instructions in terms of CPU cycles.

This package contains native Python bindings for the applying the efficient
positional population count operator to computing summary statistics for the SAM
FLAG field

## Intallation

Install with
```bash
pip3 install .
```

or locally with
```bash
python3 setup.py build_ext --inplace
```

Uninstall with
```bash
pip3 uninstall pyflagstats
```

## Example

```python
import numpy as np
import pyflagstats as fs

# Compute summary statistics for 100 million random FLAG fields.
# Completes in around 1 second.
fs.flagstats(np.random.randint(0,8192,100000000,dtype="uint16"))
```

returns (for example)

```
{'passed': array([ 624787,  312748, 2500089,  312384,  312314,  312678,  312045,
        311845, 2499502, 4999279, 2497500, 1248979,  389744,  156194,
        156029,       0], dtype=uint32), 'failed': array([ 625143,  312906, 2498840,  312818,  312129,  312802,  311869,
        312105, 2501477, 5000721, 2499178, 1249105,  390962,  155828,
        156018,       0], dtype=uint32)}
```