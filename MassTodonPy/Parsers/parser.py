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
from operator import itemgetter
from itertools import izip

def round_spec(mz, intensity, prec_digits=2):
    '''Aggregate the spectrum so that intensities of masses with the same number of significant digits are summed.'''
    mz = np.round(mz, prec_digits)
    mz, intensity = aggregate(mz, intensity)
    return mz, intensity


def trim_spectrum(mz, intensity, cut_off=100):
    '''Remove peaks below a given cut off.'''
    return mz[intensity >= cut_off], intensity[intensity >= cut_off]


def quantile_trim(mz, intensities, perc = .95):
    mz_res, intensities_res = tuple( np.array(x) for x in izip(*sorted(izip(mz, intensities), key=itemgetter(1))) )
    i = 0
    S = 0.0
    total = intensities_res.sum()
    while True:
        S += intensities_res[i]/total
        if S < 1.0-perc:
            i += 1
        else:
            break

    effective_cut_off = intensities_res[i]
    mz_res = mz_res[ intensities_res >= effective_cut_off ]
    intensities_res = intensities_res[ intensities_res >= effective_cut_off ]

    return tuple( np.array(x) for x in izip(*sorted(izip(mz_res, intensities_res), key=itemgetter(0))) ), effective_cut_off


def get_mzxml(path, prec_digits=2):
    '''Generate a sequence of rounded and trimmed spectra from individual runs of the instrument.'''
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            mz, intensity = round_spec(mz, intensity, prec_digits)
            yield mz, intensity


def read_mzxml(path, prec_digits=2):
    '''Read and merge runs of the instrument.'''
    mz, intensity = reduce( merge_runs, get_mzxml(path, prec_digits) )
    return mz, intensity


def read_txt(path, prec_digits=2):
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
    mz, intensity = round_spec(mz, intensity, prec_digits)
    return mz, intensity


def parse_path(path):
    '''Parsers path to the file.'''
    file_path, file_ext  = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext


def read_n_preprocess_spectrum( path    = None,
                                spectrum= None,
                                prec_digits = 2,
                                cut_off = None,
                                opt_P   = None  ):
    '''Read the spectrum from path or directly from provided variable, round masses to prec_digits significant digits, apply intensity based thresholding of peaks.

    It can either directly cut off peaks with intensities below the provided cut_off value or compute the smallest set of peaks with their cumulative intensity equal to opt_P percent of the total intensity.
    '''

    assert path or spectrum, "No path to mass spectrum or no mass spectrum."

    assert not (path and spectrum), "Please decide if you pass the spectrum as the argument (provide spectrum argument for method read_spectrum) or you want MassTodon to read the spectrum from file following the provided path."

    assert not (cut_off and opt_P), "Please decide if you want to apply an intensity based cut-off or to choose the optimal P-set of experimental peaks, e.g. a representative of the class of the smallest sets of peaks with a probability at least P."

    if not (cut_off or opt_P):
        print '\nAttention! \nYou did not provide a cut-off value for noise intensities, nor did you provide how probable should optimal P-set be. Thus, we will use a default value of opt_P = .99. \n\n It means that we trim the spectrum so that the remaining peaks amount for at least 99 per cent of all intensity and so that these peaks are the heighest in the spectrum.'
        opt_P = .99
    if path:
        file_path, file_name, file_ext = parse_path(path)
        file_ext = file_ext.lower()
        reader = {  '':         read_txt,
                    '.txt':     read_txt,
                    '.mzxml':   read_mzxml
        }[file_ext]
        spectrum = reader(path, prec_digits)
    else:
        spectrum = round_spec(*spectrum, prec_digits=prec_digits)

    spectra = {}
    spectra['original'] = spectrum
    spectra['original total intensity'] = sum(spectrum[1])
    if cut_off:
        spectra['trimmed'] = trim_spectrum( *spectrum, cut_off=cut_off)
    else:
        spectra['trimmed'], cut_off = quantile_trim( *spectrum, perc = opt_P)

    spectra['trimmed total intensity'] = spectra['trimmed'][1].sum()
    spectra['cut_off'] = cut_off
    return spectra
