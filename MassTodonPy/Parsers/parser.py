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
    mz_trimmed = mz[intensity <= cutOff]
    mz = mz[intensity > cutOff]
    intensity_trimmed = intensity[intensity <= cutOff]
    intensity = intensity[intensity > cutOff]
    return mz_trimmed, intensity_trimmed


def percent_trim(mz, intensity, cutOff=.999):
    '''Retrieve P percent of the highest peaks.'''
    order= np.argsort(intensity)[::-1]
    spec = np.column_stack((mz,intensity))[order,]
    totalIntensity = sum(intensity)
    spec_trimmed = spec[ np.cumsum(spec[:,1])/totalIntensity > cutOff, ]
    spec = spec[ np.cumsum(spec[:,1])/totalIntensity <= cutOff, ]
    spec.sort(0)
    spec_trimmed.sort(0)
    mz, intensity = spec[:,0], spec[:,1]
    mz_trimmed, intensity_trimmed = spec_trimmed[:,0], spec_trimmed[:,1]
    return (mz, intensity), (mz_trimmed, intensity_trimmed)


def get_mzxml(path, cutOff=100, digits=2):
    '''Generate a sequence of rounded and trimmed spectra from individual runs of the instrument.'''
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            mz, intensity = round_spec(mz, intensity, )
            yield mz, intensity


def read_mzxml(path, cutOff=100, digits=2):
    '''Read and merge runs of the instrument.'''
    mz, intensity = reduce(merge_runs, get_mzxml(path, cutOff, digits))
    return mz, intensity


def read_txt(path, cutOff=100, digits=2):
    mz = []
    intensity = []
    total_intensity = 0.0
    with open(path) as f:
        for l in f:
            l = l.split()
            I = float(l[1])
            total_intensity += I
            mz.append(float(l[0]))
            intensity.append(I)
    mz = np.array(mz)
    intensity = np.array(intensity)
    mz, intensity = round_spec(mz, intensity, digits=2)
    return mz, intensity


def parse_path(path):
    '''Parsers path to the file.'''
    file_path, file_ext  = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext


#TODO: add support for mzml files.
def readSpectrum(   path    = None,
                    spectrum= None,
                    cutOff  = 100,
                    digits  = 2,
                    P       = 1.0  ):
    if path:
        file_path, file_name, file_ext = parse_path(path)
        file_ext = file_ext.lower()
        reader = {  '':         read_txt,
                    '.txt':     read_txt,
                    '.mzxml':   read_mzxml
        }[file_ext]
        spectrum = reader(path, cutOff, digits)
    else: # spectrum provided directly
        assert spectrum != None, "what kind of non-existing spectrum is it?"

    mz, intensity = spectrum

    # print intensity.max()

    total_I = sum(intensity)
    mz_trimmed, intensity_trimmed = trim_spectrum(mz, intensity, cutOff)
    total_I_after_cut_off = sum(intensity)

    # print intensity.max()

    # if P < 1.0:
    #     (mz, intensity), (mz_perc_trimmed, intensity_perc_trimmed) = percent_trim(mz, intensity, P)
    #     mz_trimmed = np.append(mz_trimmed,mz_perc_trimmed)
    #     intensity_trimmed = np.append(intensity_trimmed, intensity_perc_trimmed)


    # print intensity.max()

    return (mz, intensity), (total_I, total_I_after_cut_off), (mz_trimmed, intensity_trimmed)
