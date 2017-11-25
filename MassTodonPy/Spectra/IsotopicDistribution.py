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

from MassTodonPy.Spectra.Measure import Measure


class IsotopicDistribution(Measure):
    """Store an isotopic distribution."""

    def __init__(self, mz=np.array([]), probability=np.array([])):
        """Initialize an isotopic distribution.

        Parameters
        ----------
        mz : numpy array
            Mass to charge ratios of the isotopic distribution.
        probability : numpy array
            Probabilities of the isotopic distribution.

        """
        self.mz = np.array(mz)
        self.probability = np.array(probability)

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def probability(self):
        """Get probabilities."""
        return self.masses

    @probability.setter
    def probability(self, probability):
        """Set probabilities."""
        self.masses = probability

    def round_masses(self, precision):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.

        """
        self.round_atoms(precision)
