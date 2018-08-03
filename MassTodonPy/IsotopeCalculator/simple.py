from __future__ import absolute_import, division, print_function
from IsoSpecPy import IsoSpecPy
from math import sqrt
import numpy as np

from MassTodonPy.Data.Constants    import infinity
from MassTodonPy.Data.get_isotopes import get_isotopic_masses, get_isotopic_probabilities



class IsotopeCalculator(object):
    def __init__(self,
                 _masses        = get_isotopic_masses(),
                 _probabilities = get_isotopic_probabilities(),
                 _isotope_DB    = {}):
        """Initialize the isotopic calculator."""
        self._masses        = _masses
        self._probabilities = _probabilities
        self._isotope_DB    = _isotope_DB
        self.mean_mass      = {}
        self.mean_variance  = {}
        for el in self._probabilities:
            self.mean_mass[el], self.mean_variance[el] = \
                get_mean_and_variance(self._masses[el],
                                      self._probabilities[el])

    def monoisotopic_mass(self,
                          formula):
        return sum(self._masses[el][0] * count
                   for el, count in Formula(formula).items())
        


def isotope_calculator(_masses        = None,
                       _probabilities = None,
                       _isotope_DB    = None):
    IC = IsotopeCalculator(_masses=None,
                           _probabilities=None,
                           _isotope_DB=None)
    return IC
