# TODO: this could be further optimized by the use of numpy
# 

from __future__ import absolute_import, division, print_function
from IsoSpecPy.IsoSpecPy import IsoSpec
from math import sqrt
import numpy as np

from MassTodonPy.Data.Constants    import infinity
from MassTodonPy.Data.get_isotopes import get_isotopic_masses, get_isotopic_probabilities
from MassTodonPy.IsotopeCalculator.Misc import get_mean_and_variance
from MassTodonPy.IsotopeCalculator.Misc import cdata2numpyarray  # TODO IsoSpec 2.0
from MassTodonPy.IsotopeCalculator.envelope import envelope



class IsotopeCalculator(object):
    def __init__(self,
                 _masses        = get_isotopic_masses(),
                 _probabilities = get_isotopic_probabilities(),
                 _isotope_DB    = {}):
        """Initialize the isotopic calculator."""
        self._masses        = _masses
        self._probabilities = _probabilities
        self._isotope_DB    = _isotope_DB
        self._mean_mass      = {}
        self._mean_variance  = {}
        for el in self._probabilities:
            self._mean_mass[el], self._mean_variance[el] = \
                get_mean_and_variance(self._masses[el],
                                      self._probabilities[el])

    def monoisotopic(self, formula, q=0, g=0):
        """Calculate monoisotopic m/z or mass of a molecule."""
        mass = sum(self._masses[el][0] * formula[el] for el in formula)
        if q > 0:
            H_mass = self._masses['H'][0]
            mz = (mass + (q + g) * H_mass) / q
            return mz
        else:
            return mass

    def mean(self, formula, q=0, g=0):
        """Calculate average m/z or mass of a molecule."""
        mass = sum(self._mean_mass[el] * formula[el] for el in formula)
        if q > 0:
            H_mass = self._mean_mass['H']
            mz = (mass + (q + g) * H_mass) / q
            return mz
        else:
            return mass

    def sd(self, formula, q, g=0):
        """Calculate standard deviation of m/z or mass for a molecule."""
        var = sum(self._mean_variance[e] * formula[e] for e in formula)
        if q > 0:
            H_var = self._mean_variance['H']
            var += H_var * (q + g)
            return sqrt(var) / q
        else:
            return sqrt(var)

    def _make_envelope(self, formula, prob, sort = True):
        el_cnt  = list(formula.values())
        el_mass = [self._masses[e]        for e in formula]
        el_prob = [self._probabilities[e] for e in formula]
        E = IsoSpecPy.IsoSpec(el_cnt, el_mass, el_prob, el_prob)
        mass, logprobability, _ = E.getConfsRaw()
        m = cdata2numpyarray(mass)
        p = np.exp(cdata2numpyarray(logprobability))
        env = envelope(m, p, sort=True)
        return env

    def __call__(self,
                 formula,
                 prob     = .999,
                 q        = 0,
                 g        = 0,
                 _memoize = False):
        if _memoize:
            env_key = (str(formula), joint_probability)
            try:
                env = self._isotope_DB[env_key]
            except KeyError:
                env = self._make_envelope(formula, prob)
                self._isotope_DB[env_key] = env
        else:
            env = self._make_envelope(formula, prob)




def isotope_calculator(_masses        = get_isotopic_masses(),
                       _probabilities = get_isotopic_probabilities(),
                       _isotope_DB    = {}):
    IC = IsotopeCalculator(_masses        = _masses,
                           _probabilities = _probabilities,
                           _isotope_DB    = _isotope_DB)
    return IC
