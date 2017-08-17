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
import os

from pyteomics import mzxml # >= 3.41

from MassTodonPy.Spectra.operations import aggregate, merge_runs, trim_spectrum, round_spectrum, remove_lower_quantile


def get_mzxml(path, prec_digits=2):
    """Generate a sequence of rounded and trimmed spectra from individual runs of the instrument.

    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    prec_digits : float
        The number of digits after which the floats get rounded.

    Returns
    -------
    out : generator
        Generates tuples of numpy arrays corresponding to different runs of the experimental spectrum.
    """
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            mz, intensity = round_spectrum(mz, intensity, prec_digits)
            yield mz, intensity


def read_mzxml(path, prec_digits=2):
    """Read spectrum form an mzXml and merge runs of the instrument.

    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    prec_digits : float
        The number of digits after which the floats get rounded.

    Returns
    -------
    out : generator
        Generates tuples of numpy arrays corresponding to different runs of the experimental spectrum.
    """
    mz, intensity = reduce( merge_runs, get_mzxml(path, prec_digits) )
    return mz, intensity


def read_txt(path, prec_digits=2):
    """Read spectrum from a text file.

    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    prec_digits : float
        The number of digits after which the floats get rounded.

    Returns
    -------
    out : generator
        Generates tuples of numpy arrays corresponding to different runs of the experimental spectrum.
    """
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
    mz, intensity = round_spectrum(mz, intensity, prec_digits)
    return mz, intensity


def parse_path(path):
    """Parse path to the file.

    Parameters
    ----------
    path : str
        Any path.

    Returns
    -------
    out : tuple
        Path of file, name of file, and file's extension.
    """
    file_path, file_ext  = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext


def read_n_preprocess_spectrum( path    = None,
                                spectrum= None,
                                prec_digits = 2,
                                cut_off = None,
                                opt_P   = None  ):
    """Read the spectrum, round it, apply intensity based thresholding of peaks.

    Parameters
    ----------
    path : str
        Path to spectrum.
    spectrum : tuple
        The experimental spectrum, a tuple of m/z numpy array and intensities numpy array.
    prec_digits : float
        The number of digits after which the floats get rounded.
    cut_off : float
        The cut off value for peak intensity.
    opt_P :
        The percentage of the heighest peaks being used in the analysis.

    Warning
    -------
    You should provide either the path or the spectrum, not both simultaneously.

    Returns
    -------
    spectra : dict
        A dictionary containing the experimental spectrum split into several subspectra.

    Notes
    -----
    **spectra** contains
        * the original experimental spectrum
        * its total intensity
        * spectrum after trimming, *trimmed*
        * total intensity after trim
        * the effective peak height cut off
    """

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
        spectrum = round_spectrum( *spectrum, prec_digits=prec_digits)

    spectra = {}
    spectra['original'] = spectrum
    spectra['original total intensity'] = sum(spectrum[1])

    if opt_P: # cut_off == None
        cut_off = remove_lower_quantile(*spectrum, retained_percentage = opt_P)

    spectra['trimmed'], spectra['left out'] = trim_spectrum( *spectrum, cut_off=cut_off)
    spectra['total intensity after trim'] = spectra['trimmed'][1].sum()
    spectra['trimmed intensity'] = spectra['original total intensity'] - spectra['total intensity after trim']
    spectra['cut_off'] = cut_off

    return spectra
