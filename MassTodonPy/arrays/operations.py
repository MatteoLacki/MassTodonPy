"""Various operations on arrays."""

import numpy as np

def is_sorted(x):
    return np.all(np.diff(x))
