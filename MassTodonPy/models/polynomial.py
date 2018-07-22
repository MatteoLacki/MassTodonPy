"""Class that wraps up the polynomial fitting."""

import numpy as np

from MassTodonPy.models.two_dimensional import Model2D


class Polynomial(Model2D):
    def __init__(self, degree = 3):
        """Initialize the Polynomial class.

        Parameters
        ----------
        degree : int
            The degree of the fitted polynomial.
        """
        assert degree >= 1, "The degree must be greater or equal to 1."
        self.degree = int(3)

    def fit(self, x, y):
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
        self.x, self.y   = x, y
        self.x_min = min(self.x)
        self.x_max = max(self.x)
        self._coefs      = np.polyfit(self.x, self.y, self.degree)
        self._polynomial = np.poly1d(self._coefs)

    def coef(self):
        """Get the coefficients of the polynomial fitted with least squares."""
        return self._coefs

    def __call__(self, x):
        """Extrapolate the value of the fitted polynomial at 'x'.

        Parameters
        ----------
        x : np.array
            Values at which to extrapolate the polynomial.

        """
        return self._polynomial(x)

    def __repr__(self):
        return "Ich bin ein Polynom. Hast du Angst?"


def polynomial(x, y, degree = 3):
    """Fit a scipy-spline to the data.

    Parameters:
    x : np.array
        1D control.
    y : np.array
        1D response.
    degree : int
        The degree of the fitted polynomial.
    """
    p = Polynomial(degree)
    p.fit(x, y)
    return p