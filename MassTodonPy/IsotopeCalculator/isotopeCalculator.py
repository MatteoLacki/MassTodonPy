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
import numpy as np
from collections import Counter
from time import time
from math import sqrt
from IsoSpecPy import IsoSpecPy

from MassTodonPy.Data.get_isotopes import get_isotopic_masses_and_probabilities
from MassTodonPy.Parsers.formula_parser import parse_formula
from MassTodonPy.Spectra.Spectrum import TheoreticalSpectrum
from MassTodonPy.Spectra.operations import cdata2numpyarray,\
                                           aggregate,\
                                           aggregate_envelopes


def get_mean_and_variance(X, weights):
    X = np.array(X)
    probs = np.array(weights)
    probs = probs/sum(probs)
    average = np.dot(X, probs)
    variance = np.dot((X-average)**2, probs)
    return average, variance


def check_charges(q, g):
    assert isinstance(q, int) and q > 0, "q must be a positive integer."
    assert isinstance(g, int) and g >= 0, "g must be a non-negative integer."


# TODO: convolute spectra with diffs spectra instead of Dirac deltas.
class IsotopeCalculator:
    """A class for isotope calculations."""

    iso_masses, iso_probs = get_isotopic_masses_and_probabilities()

    def __init__(self,
                 prec_digits=2,
                 iso_masses=None,
                 iso_probs=None,
                 verbose=False):
        """Instantiate the isotopic calculator."""

        if iso_masses is not None and iso_probs is not None:
            self.iso_masses = iso_masses
            self.iso_probs = iso_probs
        self.mean_mass = {}
        self.mean_variance = {}
        for el in self.iso_probs:
            self.mean_mass[el], self.mean_variance[el] = \
                get_mean_and_variance(self.iso_masses[el],
                                      self.iso_probs[el])
        self.isotope_DB = {}  # TODO replace with the IsoSpec generator.
        self.prec_digits = prec_digits
        self.verbose = verbose
        self.stats = Counter()  # TODO replace with logger

    def get_monoisotopic_mz(self, formula, q, g=0):
        """Calculate monoisotopic mass of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        mass = sum(self.iso_masses[el][0]*elCnt
                   for el, elCnt in atomCnt.items())
        hydrogen_mass = self.iso_masses['H'][0]
        mz = (mass + (q + g)*hydrogen_mass) / q
        return mz

    def get_mean_mz(self, formula, q, g=0):
        """Calculate average mass of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        mass = sum(self.mean_mass[el]*elCnt for el, elCnt in atomCnt.items())
        hydrogen_mass = self.mean_mass['H']
        mz = (mass + (q + g)*hydrogen_mass)/q
        return mz

    def get_mz_sd(self, formula, q, g=0):
        """Calculate m/z standard deviation of a molecule."""
        check_charges(q, g)
        atomCnt = parse_formula(formula)
        sd = sqrt(sum(self.mean_variance[el]*elCnt
                      for el, elCnt in atomCnt.items()))
        return sd/q

    def __make_envelope(self, formula, joint_probability, memoize=False):
        T0 = time()
        counts = []
        isotope_masses = []
        isotope_probs = []
        atomCnt = parse_formula(formula)

        for el, cnt in atomCnt.items():
            counts.append(cnt)
            isotope_masses.append(self.iso_masses[el])
            isotope_probs.append(self.iso_probs[el])

        envelope = IsoSpecPy.IsoSpec(counts,
                                     isotope_masses,
                                     isotope_probs,
                                     joint_probability)

        masses, logprobs, _ = envelope.getConfsRaw()
        masses = cdata2numpyarray(masses)
        probs = np.exp(cdata2numpyarray(logprobs))
        masses, probs = aggregate_envelopes(masses, probs, self.prec_digits)

        if memoize:
            # TODO get rid of it when IsoSpec will be fastest possible.
            self.isotope_DB[(formula, joint_probability)] = (masses, probs)

        T1 = time()
        self.stats['Envelopes Generation Total T'] += T1-T0
        return masses.copy(), probs.copy()

    def get_envelope(self,
                     joint_probability,
                     formula,
                     q=0,
                     g=0,
                     memoize=False):
        """Get an isotopic envelope consisting of a numpy array
        of masses and numpy array of probabilities.

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

        hydrogen_mass = self.mean_mass['H']
        masses = np.around((masses + (g + q)*hydrogen_mass)/q,
                           decimals=self.prec_digits)

        masses, probs = aggregate(masses, probs)
        return TheoreticalSpectrum(mass=masses, probability=probs)
