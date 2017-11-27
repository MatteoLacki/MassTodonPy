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
from __future__ import absolute_import, division, print_function
from collections import defaultdict
from IsoSpecPy import IsoSpecPy
from math import fsum
from math import sqrt
import numpy as np

from MassTodonPy.Data.get_isotopes import get_isotopic_masses_and_probabilities
from MassTodonPy.Formula.Formula import Formula
from MassTodonPy.IsotopeCalculator.Misc import cdata2numpyarray
from MassTodonPy.IsotopeCalculator.Misc import get_mean_and_variance
from MassTodonPy.IsotopeCalculator.Misc import check_charges
from MassTodonPy.MoleculeMaker.formula_parser import get_formula_parser
from MassTodonPy.Spectra.IsotopicDistribution import IsotopicDistribution

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


# convolute spectra with diffs spectra instead of Dirac deltas.
class IsotopeCalculator(object):
    """Perform isotope calculations.

    Parameters
    ==========
    mz_precision : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    joint_probability : float
        The joint probability threshold for generating
        theoretical isotopic distributions.
    _isotope_masses : dict
        The isotopic masses, e.g. from IUPAC.
    _isotope_probabilities : dict
        The isotopic frequencies, e.g. from IUPAC.

    """
    _isotope_masses, _isotope_probabilities =\
        get_isotopic_masses_and_probabilities()
    _isotope_DB = {}  # replace with the IsoSpec generator.

    def __init__(self,
                 mz_precision=2,
                 joint_probability=.999,
                 _isotope_masses=None,
                 _isotope_probabilities=None,
                 _isotope_DB=None):
        """Initialize the isotopic calculator."""
        if _isotope_masses:
            self._isotope_masses = _isotope_masses
        if _isotope_probabilities:
            self._isotope_probabilities = _isotope_probabilities
        if _isotope_DB:
            self._isotope_DB = _isotope_DB
        self.mean_mass = {}
        self.mean_variance = {}
        for el in self._isotope_probabilities:
            self.mean_mass[el], self.mean_variance[el] = \
                get_mean_and_variance(self._isotope_masses[el],
                                      self._isotope_probabilities[el])
        self.mz_precision = mz_precision
        self.joint_probability = joint_probability

    def get_monoisotopic_mz(self, formula, q, g=0):
        """Calculate monoisotopic mass of a molecule."""
        check_charges(q, g)
        mass = sum(self._isotope_masses[el][0] * count
                   for el, count in Formula(formula).items())
        hydrogen_mass = self._isotope_masses['H'][0]
        mz = (mass + (q + g) * hydrogen_mass) / q
        return mz

    def get_mean_mz(self, formula, q, g=0):
        """Calculate average mass of a molecule."""
        check_charges(q, g)
        mass = sum(self.mean_mass[el] * count
                   for el, count in Formula(formula).items())
        hydrogen_mass = self.mean_mass['H']
        mz = (mass + (q + g) * hydrogen_mass) / q
        return mz

    def get_mz_sd(self, formula, q, g=0):
        """Calculate m/z standard deviation of a molecule."""
        check_charges(q, g)
        sd = sqrt(sum(self.mean_variance[el] * elCnt
                      for el, elCnt in Formula(formula).items()))
        return sd / q

    def __make_envelope(self,
                        formula,
                        joint_probability,
                        memoize=False):
        counts = []
        isotope_masses = []
        isotope_probs = []
        for el, cnt in Formula(formula).items():
            counts.append(cnt)
            isotope_masses.append(self._isotope_masses[el])
            isotope_probs.append(self._isotope_probabilities[el])
        envelope = IsoSpecPy.IsoSpec(counts,
                                     isotope_masses,
                                     isotope_probs,
                                     joint_probability)
        masses, logprobs, _ = envelope.getConfsRaw()
        masses = cdata2numpyarray(masses)  # TODO IsoSpec 2.0
        probs = np.exp(cdata2numpyarray(logprobs))
        # TODO : get rid of this when IsoSpec 2.0 is in place
        masses, probs = aggregate_envelopes(masses,
                                            probs,
                                            self.mz_precision)
        if memoize:
            # TODO get rid of it when IsoSpec will be fastest possible.
            self._isotope_DB[(str(formula), joint_probability)] = (masses, probs)
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
        formula : str or Formula object
            E.g. 'C63H98N18O13S1', 'H2O'
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
            masses, probs = self._isotope_DB[(str(formula),
                                              joint_probability)]
            masses, probs = masses.copy(), probs.copy()
        except KeyError:
            masses, probs = self.__make_envelope(formula,
                                                 joint_probability,
                                                 memoize)

        if joint_probability is None:
            joint_probability = self.joint_probability

        hydrogen_mass = self.mean_mass['H']
        masses = np.around((masses + (g + q) * hydrogen_mass) / q,
                           decimals=self.mz_precision)

        masses, probs = aggregate(masses, probs)
        return IsotopicDistribution(mz=masses, probability=probs)
