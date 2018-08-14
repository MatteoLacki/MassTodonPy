import numpy as np


def by_rounding(x, y, digits):
    x_r    = np.around(x, digits)
    x_r, i = np.unique(x_r, return_inverse=True)
    y_agg  = np.bincount(i, weights=y)
    return x_r, y_agg
