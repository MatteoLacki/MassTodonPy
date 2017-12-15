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

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Spectra.Measure import Measure


class ExperimentalSpectrum(Measure):
    """Store an experimental spectrum."""

    def __init__(self, mz=np.array([]), intensity=np.array([])):
        """Initialize an experimental spectrum.

        Parameters
        ----------
        mz : numpy array
            Spectral mass to charge ratios.
        intensity : numpy array
            Spectral intensities.

        """
        self.mz = np.array(mz)
        self.intensity = np.array(intensity)
        self._store_names = ('m/z', 'intensity')

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def intensity(self):
        """Get intensities."""
        return self.masses

    @intensity.setter
    def intensity(self, intensity):
        """Set intensities."""
        self.masses = intensity

    def round_mz(self, precision=infinity):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
            Defaults to 'inf', which prevents any rounding.

        """
        self.round_atoms(precision)

    def total_intensity(self):
        return self.intensity.sum()

    def trim_intensity(self, cut_off):
        """Trim intensities below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        self.trim(cut_off)
