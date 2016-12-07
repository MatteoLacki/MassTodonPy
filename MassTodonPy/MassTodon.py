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

from Formulator         import makeFragments
from IsotopeCalculator  import isotopeCalculator
from PeakPicker         import peakPicker

class MassTodon():
    def __init__( self, fasta, precursorCharge, modifications={}, massPrecDigits = 3, isoMasses=None, isoProbs=None ):
        self.fasta  = fasta
        self.Q      = precursorCharge
        self.isoCalc= isotopeCalculator(massPrecDigits, isoMasses, isoProbs)

    def randomSpectrum(
            self,
            ionsNo,
            fragScheme      = 'cz',
            aaPerOneCharge  = 5,
            jointProb       = .999,
            scale           = .01,
            percentPeaks    =.2     ):

        masses, intensities = self.isoCalc.randomSpectrum(
            self.fasta,
            self.Q,
            ionsNo,
            fragScheme='cz',
            aaPerOneCharge=5,
            jointProb=.999,
            scale =.01          )

        noise_masses, noise_intensities = self.isoCalc.addNoise(
            masses,
            intensities,
            percentPeaks        )

        return masses, intensities, noise_masses, noise_intensities

    def setMassSpectrum(self, massSpectrum):
        self.massSpectrum = massSpectrum

    def pickPeaks(self):
        self.peakPicker = PeakPicker()
