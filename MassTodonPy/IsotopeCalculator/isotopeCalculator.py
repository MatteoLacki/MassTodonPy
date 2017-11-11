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

import numpy as np

from IsoSpecPy import IsoSpecPy
from collections import Counter
from time import time

from MassTodonPy.Data.get_data import get_isotopic_masses_and_probabilities
from MassTodonPy.Parsers.formula_parser import formulaParser
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

    def __init__(self,
                 prec_digits=2,
                 iso_masses=None,
                 iso_probs=None,
                 verbose=False):
        """Instantiate the isotopic calculator."""

        if iso_masses is None or iso_probs is None:
            self.iso_masses, self.iso_probs = \
                get_isotopic_masses_and_probabilities()

        self.iso_means_vars = {
            el: get_mean_and_variance(self.iso_masses[el],
                                      self.iso_probs[el])
            for el in self.iso_probs}

        self.isotope_DB = {}
        self.prec_digits = prec_digits
        self.formParser = formulaParser()
        self.verbose = verbose
        self.stats = Counter()

    def __getMonoisotopicMass(self, atomCnt):
        """Calculate monoisotopic mass of an atom count."""
        return sum(self.iso_masses[el][0]*elCnt
                   for el, elCnt in atomCnt.items())  # .items is pythonic!

    def __getMassMean(self, atomCnt):
        """Calculate average mass of an atom count."""
        return sum(self.iso_means_vars[el][0]*elCnt
                   for el, elCnt in atomCnt.items())

    def __getMassVar(self, atomCnt):
        """Calculate mass variance of an atom count."""
        return sum(self.iso_means_vars[el][1]*elCnt
                   for el, elCnt in atomCnt)

    def __getSummary(self, atomCnt_str):
        atomCnt = self.formParser.parse(atomCnt_str)
        return (self._getMonoisotopicMass(atomCnt),
                self._getMassMean(atomCnt),
                self._getMassVar(atomCnt))

    def __getOldEnvelope(self, atomCnt_str, jP, prec_digits):
        masses, probs = self.isotope_DB[(atomCnt_str, jP, prec_digits)]
        return masses.copy(), probs.copy()

    def __make_envelope(self, formula, joint_probability):
        T0 = time()
        counts = []
        isotope_masses = []
        isotope_probs = []
        atomCnt = self.formParser.parse(formula.formula)

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
        self.isotope_DB[(formula.formula, joint_probability)] = (masses, probs)
        T1 = time()

        self.stats['Envelopes Generation Total T'] += T1-T0
        return masses.copy(), probs.copy()

    def get_envelope(self,
                     formula,
                     joint_probability):
        """Get an isotopic envelope consisting of a numpy array
        of masses and numpy array of probabilities.

        Parameters
        ----------
        atomCnt_str : str
            The chemical formula of a molecular species.
        formula : a namedtuple
            e.g. Formula(name='precursor', formula='C63H98N18O13S1', q=3, g=0)
        joint_probability : float
            The joint probability of the theoretical isotopic envelope.

        Returns
        -------
        out : tuple
            A tuple containing the theoretical spectrum:
            mass over charge values and intensities, both numpy arrays.
        """
        try:
            masses, probs = self.isotope_DB[(formula.formula,
                                             joint_probability)]

            masses, probs = masses.copy(), probs.copy()
        except KeyError:
            masses, probs = self.__make_envelope(formula,
                                                 joint_probability)

        if formula.q is not 0:
            # TODO get proper mass of a hydrogen!
            masses = np.around(
                (masses + formula.g + formula.q)/formula.q,
                decimals=self.prec_digits)
        masses, probs = aggregate(masses, probs)
        return masses, probs
