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

        self.isotope_DB = {}
        self.prec_digits = prec_digits
        self.verbose = verbose
        self.stats = Counter()

    def get_monoisotopic_mz(self, molecule):
        """Calculate monoisotopic mass of a molecule."""
        atomCnt = parse_formula(molecule.formula)
        mass = sum(self.iso_masses[el][0]*elCnt
                   for el, elCnt in atomCnt.items())
        hydrogen_mass = self.iso_masses['H'][0]
        mz = (mass + (molecule.q + molecule.g)*hydrogen_mass)/molecule.q
        return mz

    def get_mean_mz(self, molecule):
        """Calculate average mass of a molecule."""
        atomCnt = parse_formula(molecule.formula)
        mass = sum(self.mean_mass[el]*elCnt for el, elCnt in atomCnt.items())
        hydrogen_mass = self.mean_mass['H']
        mz = (mass + (molecule.q + molecule.g)*hydrogen_mass)/molecule.q
        return mz

    def get_mz_sd(self, molecule):
        """Calculate m/z standard deviation of a molecule."""
        atomCnt = parse_formula(molecule.formula)
        sd = sqrt(sum(self.mean_variance[el]*elCnt
                      for el, elCnt in atomCnt.items()))
        return sd/molecule.q

    def __make_envelope(self, formula, joint_probability):
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

        # memoization
        # TODO get rid of it when IsoSpec using stats convolution
        self.isotope_DB[(formula, joint_probability)] = (masses, probs)
        T1 = time()

        self.stats['Envelopes Generation Total T'] += T1-T0
        return masses.copy(), probs.copy()

    def get_envelope(self,
                     molecule,
                     joint_probability):
        """Get an isotopic envelope consisting of a numpy array
        of masses and numpy array of probabilities.

        Parameters
        ----------
        atomCnt_str : str
            The chemical formula of a molecular species.
        molecule : a namedtuple
            e.g. Molecule(name='precursor', formula='C63H98N18O13S1', q=3, g=0)
        joint_probability : float
            The joint probability of the theoretical isotopic envelope.

        Returns
        -------
        out : tuple
            A tuple containing the theoretical spectrum:
            mass over charge values and intensities, both numpy arrays.
        """
        try:
            masses, probs = self.isotope_DB[(molecule.formula,
                                             joint_probability)]

            masses, probs = masses.copy(), probs.copy()
        except KeyError:
            masses, probs = self.__make_envelope(molecule.formula,
                                                 joint_probability)

        if molecule.q is not 0:
            # TODO get proper mass of a hydrogen!
            hydrogen_mass = self.mean_mass['H']
            masses = np.around(
                (masses + (molecule.g + molecule.q)*hydrogen_mass)/molecule.q,
                decimals=self.prec_digits)
        masses, probs = aggregate(masses, probs)
        return masses, probs
