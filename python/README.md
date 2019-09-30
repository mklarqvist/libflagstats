# pyflagstats

python bindings for libflagstats. 

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

fs.flagstats(np.random.randint(0,8192,10000000,dtype="uint16"))
```

returns (for example)

```
{'passed': array([ 624787,  312748, 2500089,  312384,  312314,  312678,  312045,
        311845, 2499502, 4999279, 2497500, 1248979,  389744,  156194,
        156029,       0], dtype=uint32), 'failed': array([ 625143,  312906, 2498840,  312818,  312129,  312802,  311869,
        312105, 2501477, 5000721, 2499178, 1249105,  390962,  155828,
        156018,       0], dtype=uint32)}
```