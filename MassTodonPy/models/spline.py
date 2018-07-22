import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LSQUnivariateSpline

from MassTodonPy.models.two_dimensional import Model2D
from MassTodonPy.arrays.operations import dedup_sort

class Spline(Model2D):
    def fit(self, x, y,
            drop_duplicates = True,
            sort            = True):
        """Fit a spline to 2D data.

        Parameters
        ----------
        x : np.array
            One-dimensional control variables.
        y : np.array
            One-dimensional response variable.
        drop_duplicates : logical
            Drop duplicate rows in 'xy' data-frame.
        sort : logical
            Sort the 'xy' data-frame by 'x'.
        Returns
        -------
        tuple of np.arrays: fixed 'x' and 'y'.
        """
        self.x, self.y = dedup_sort(x, y, drop_duplicates, sort)
        t = self.x_percentiles(len(self.x)//1000)
        self.x_min = t[0]
        self.x_max = t[-1]
        t = t[1:-1]
        self._spline = LSQUnivariateSpline(self.x, self.y, t)

    def __call__(self, x):
        return self._spline(x)

    def __repr__(self):
        return "Ich bin ein Spline. Hast du Angst?"

def spline(x, y, 
           drop_duplicates = True,
           sort            = True):
    """Fit a scipy-spline to the data.

    Arguments:
        x : np.array
        y : np.array
    """
    s = Spline()
    s.fit(x, y, drop_duplicates, sort)
    return s