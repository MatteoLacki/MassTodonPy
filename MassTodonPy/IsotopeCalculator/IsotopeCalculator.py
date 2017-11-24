# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.
from __future__ import division
from collections import Counter
from collections import defaultdict
from IsoSpecPy import IsoSpecPy
from math import fsum
from math import sqrt
import numpy as np
from six.moves import range
from time import time

from MassTodonPy.Data.get_isotopes import get_isotopic_masses_and_probabilities
from MassTodonPy.MoleculeMaker.formula_parser import parse_formula
from MassTodonPy.Spectra.Operations import aggregate
from MassTodonPy.Spectra.Spectra import TheoreticalSpectrum as TheoSpec


def cdata2numpyarray(x):
    """Turn c-data into a numpy array.

    Parameters
    ----------
    x : cdata table
        A table of cdata from cffi.
    Returns
    -------
    res : array
        A numpy array of numbers.
    """
    res = np.empty(len(x))
    for i in range(len(x)):
        res[i] = x[i]
    return res


# This is used by the IsotopeCalculator.
# get rid of this when IsoSpec 2.0 is in place
def aggregate_envelopes(masses, probs, digits=2):
    """Aggregate theoretical envelopes.

    Parameters
    ----------
    masses : array
        An array of isotopologues' masses.
    probs : array
        An array of isotopologues' probabilities.
    digits : int
        The number of significant digits used
        while rounding the masses of isotopologues.

    Returns
    ----------
    out : tuple
        A theoretical spectrum of a given resolution.
    """
    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array([k for k in lists])
    newProbs = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs, op_flags=['readwrite']),
                          newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs


def get_mean_and_variance(X, weights):
    """Get mean and variance of X."""
    X = np.array(X)
    probs = np.array(weights)
    probs = probs / sum(probs)
    average = np.dot(X, probs)
    variance = np.dot((X - average) ** 2, probs)
    return average, variance


def check_charges(q, g):
    """Assert q and g are integers, q is positive, g nonnegative."""
    assert isinstance(q, int) and q > 0, "q must be a positive integer."
    assert isinstance(g, int) and g >= 0, "g must be a non-negative integer."


# convolute spectra with diffs spectra instead of Dirac deltas.
class IsotopeCalculator(object):
    """Perform isotope calculations.

    Parameters
    ==========
    _isotope_masses : dict
        The isotopic masses, e.g. from IUPAC.

    _isotope_probabilities : dict
        The isotopic frequencies, e.g. from IUPAC.
    """

    _isotope_masses, _isotope_probabilities =\
        get_isotopic_masses_and_probabilities()

    def __init__(self,
                 mz_precision_digits=2,
                 joint_probability=.999,
                 _isotope_masses=None,
                 _isotope_probabilities=None,
                 _verbose=False):
        """Initialize the isotopic calculator."""
        if _isotope_masses is not None:
            self._isotope_masses = _isotope_masses
        if _isotope_probabilities is not None:
            self._isotope_probabilities = _isotope_probabilities
        self.mean_mass = {}
        self.mean_variance = {}
        for el in self._isotope_probabilities:
            self.mean_mass[el], self.mean_variance[el] = \
                get_mean_and_variance(self._isotope_masses[el],
                                      self._isotope_probabilities[el])
        self.isotope_DB = {}  # replace with the IsoSpec generator.
        self.mz_precision_digits = mz_precision_digits
        self._verbose = _verbose
        self.joint_probability = joint_probability
        self.stats = Counter()  # logger

    def get_monoisotopic_mz(self, formula, q, g=0):
        """Calculate monoisotopic mass of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        mass = sum(self._isotope_masses[el][0] * elCnt
                   for el, elCnt in atomCnt.items())
        hydrogen_mass = self._isotope_masses['H'][0]
        mz = (mass + (q + g) * hydrogen_mass) / q
        return mz

    def get_mean_mz(self, formula, q, g=0):
        """Calculate average mass of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        mass = sum(self.mean_mass[el] * elCnt for el, elCnt in atomCnt.items())
        hydrogen_mass = self.mean_mass['H']
        mz = (mass + (q + g) * hydrogen_mass) / q
        return mz

    def get_mz_sd(self, formula, q, g=0):
        """Calculate m/z standard deviation of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        sd = sqrt(sum(self.mean_variance[el] * elCnt
                      for el, elCnt in atomCnt.items()))
        return sd / q

    def __make_envelope(self, formula, joint_probability, memoize=False):
        T0 = time()
        counts = []
        isotope_masses = []
        isotope_probs = []
        atomCnt = parse_formula(formula)

        for el, cnt in atomCnt.items():
            counts.append(cnt)
            isotope_masses.append(self._isotope_masses[el])
            isotope_probs.append(self._isotope_probabilities[el])

        envelope = IsoSpecPy.IsoSpec(counts,
                                     isotope_masses,
                                     isotope_probs,
                                     joint_probability)

        masses, logprobs, _ = envelope.getConfsRaw()
        masses = cdata2numpyarray(masses)
        probs = np.exp(cdata2numpyarray(logprobs))

        # TODO : get rid of this when IsoSpec 2.0 is in place
        masses, probs = aggregate_envelopes(masses,
                                            probs,
                                            self.mz_precision_digits)

        if memoize:
            # TODO get rid of it when IsoSpec will be fastest possible.
            self.isotope_DB[(formula, joint_probability)] = (masses, probs)

        T1 = time()
        self.stats['Envelopes Generation Total T'] += T1 - T0
        return masses.copy(), probs.copy()

    def get_envelope(self,
                     formula,
                     joint_probability=None,
                     q=0,
                     g=0,
                     memoize=False):
        """Get an isotopic envelope.

        Parameters
        ----------
        formula : str
            E.g. 'C63H98N18O13S1'
        joint_probability : float
            The joint probability of the theoretical isotopic envelope.
        q : int
            The charge state of the molecule.
        g : int
            Additional hydrogen atoms:
            referred to as quenched charge in the MassTodon paper.
        memoize : bool
            Should the call be memoized: by default no.
        Returns
        -------
        out : tuple
            A tuple containing the theoretical spectrum:
            mass over charge values and intensities, both numpy arrays.
        """
        check_charges(q, g)

        try:
            masses, probs = self.isotope_DB[(formula, joint_probability)]
            masses, probs = masses.copy(), probs.copy()
        except KeyError:
            masses, probs = self.__make_envelope(formula,
                                                 joint_probability,
                                                 memoize)

        if joint_probability is None:
            joint_probability = self.joint_probability

        hydrogen_mass = self.mean_mass['H']
        masses = np.around((masses + (g + q) * hydrogen_mass) / q,
                           decimals=self.mz_precision_digits)

        masses, probs = aggregate(masses, probs)
        return TheoSpec(mz=masses, probability=probs)
