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

from pyteomics          import mzxml # >= 3.41
from IsotopeCalculator  import aggregate
import numpy as np
import os

def merge_runs(spec1, spec2):
    mz = np.concatenate((spec1[0], spec2[0]))
    I  = np.concatenate((spec1[1], spec2[1]))
    return aggregate(mz, I)

def get_spectra(path, cutOff=100, digits=2):
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz          = spectrum['m/z array']
            intensity   = spectrum['intensity array']
            mz          = mz[intensity > cutOff]
            intensity   = intensity[intensity > cutOff]
            mz          = np.round(mz, digits)
            mz, intensity= aggregate(mz, intensity)
            yield (mz, intensity)

def readMzXml(path, cutOff=100, digits=2):
    mz, intensity = reduce(merge_runs, get_spectra(path, cutOff, digits))
    return mz, intensity

def readSpectrum(path, cutOff=100, digits=2):
    file_path, ext = os.path.splitext(path)
    ext = ext.lower()
    reader = {
        '.mzxml':readMzXml
    }[ext]
    mz, intensity = reader(path, cutOff, digits)
    return mz, intensity
