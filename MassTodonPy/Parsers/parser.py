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

from pyteomics   import mzxml # >= 3.41
from collections import defaultdict
from MassTodonPy.IsotopeCalculator import aggregate, merge_runs
import numpy as np
import os

def round_spec(mz, intensity, digits=2):
    '''Aggregate the spectrum so that intensities of masses with the same number of significant digits are summed.'''
    mz = np.round(mz, digits)
    mz, intensity = aggregate(mz, intensity)
    return mz, intensity


def trim_spectrum(mz, intensity, cutOff=100):
    '''Remove peaks below a given cut off.'''
    mz = mz[intensity > cutOff]
    intensity = intensity[intensity > cutOff]
    return mz, intensity


def percent_trim(mz, intensity, cutOff=.999):
    '''Retrieve P percent of the highest peaks.'''
    order= np.argsort(intensity)[::-1]
    spec = np.column_stack((mz,intensity))[order,]
    totalIntensity = sum(intensity)
    spec = spec[ np.cumsum(spec[:,1])/totalIntensity < cutOff, ]
    spec.sort(0)
    mz, intensity = spec[:,0], spec[:,1]
    return mz, intensity


def get_mzxml(path, cutOff=100, digits=2):
    '''Generate a sequence of rounded and trimmed spectra from individual runs of the instrument.'''
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            if cutOff:
                mz, intensity = trim_spectrum(mz, intensity, cutOff)
            mz, intensity = round_spec(mz, intensity, )
            yield (mz, intensity)


def read_mzxml(path, cutOff=100, digits=2):
    '''Read and merge runs of the instrument.'''
    mz, intensity = reduce(merge_runs, get_mzxml(path, cutOff, digits))
    return mz, intensity


def read_txt(path, cutOff=100, digits=2):
    mz = []
    intensity = []
    with open(path) as f:
        for l in f:
            l = l.split()
            I = float(l[1])
            if I>= cutOff:
                mz.append(float(l[0]))
                intensity.append(I)
    mz = np.array(mz)
    intensity = np.array(intensity)
    mz, intensity = round_spec(mz, intensity, digits=2)
    return mz, intensity

#TODO: add support for mzml files.
def readSpectrum(path, cutOff=100, digits=2, P=1.0):
    file_path, file_ext  = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1])

    file_ext = file_ext.lower()
    reader = {  '':         read_txt,
                '.txt':     read_txt,
                '.mzxml':   read_mzxml
    }[file_ext]
    mz, intensity = reader(path, cutOff, digits)
    if P < 1.0:
        mz, intensity = percent_trim(mz, intensity, P)
    return file_path, file_name, (mz, intensity)
