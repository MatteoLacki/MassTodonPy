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

# from MassTodonPy.Spectra.operations import trim_spectrum,\
#                                            remove_lower_quantile
import os
from MassTodonPy.Spectra.Spectra import Spectrum,\
    read_spectrum_from_txt,\
    read_spectrum_from_mzxml,\
    threshold_n_round_n_aggregate


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
    file_path, file_ext = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext


def read_n_preprocess_spectrum(
        spectrum='',
        spectral_intensity_cut_off=0.0,
        percentage_of_heighest_peaks_used=1.0,
        precision_digits=2,
        _verbose=False):
    """Read the spectrum, round it, apply intensity based thresholding of peaks.

    Parameters
    ----------
    spectrum : path string or tuple
        Either the path to the spectrum or the experimental spectrum itself,
        in form of a tuple of m/z numpy array and intensities numpy array.
        The path can end up with extension:
            *.txt, for spectra saved in tab separated format.
            *.mzXml, for spectra saved with mxXml format.

    spectral_intensity_cut_off : float
        The cut off value for peak intensity.

    percentage_of_heighest_peaks_used : float
        The percentage of the heighest peaks being used in the analysis.

    precision_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.

    Returns
    -------
    spectra : dict
        A dictionary containing the experimental spectrum
        split into several subspectra.

    Notes
    -----
    **spectra** contains
        * the original experimental spectrum
        * its total intensity
        * spectrum after trimming, *trimmed*
        * total intensity after trim
        * the effective peak height cut off
    """

    assert spectrum is not "", "Provide a spectrum to analyze!"

    assert isinstance(spectrum, (tuple, str, Spectrum)),\
        "spectrum is not a string nor a tuple."

    assert isinstance(spectral_intensity_cut_off, (int, float)) and\
        spectral_intensity_cut_off >= 0,\
        'Intensity cut off should be nonnegative.'

    assert isinstance(percentage_of_heighest_peaks_used, float) and\
        0.0 <= percentage_of_heighest_peaks_used <= 1.0,\
        'Percentage cut off should take values between 0 and 1.'

    assert isinstance(precision_digits, int) and precision_digits >= 0,\
        'Precision digits should be a natural number or zero.'

    if isinstance(spectrum, str):
        file_path, file_name, file_ext = parse_path(spectrum)
        file_ext = file_ext.lower()
        try:
            reader = {'': read_spectrum_from_txt,
                      '.txt': read_spectrum_from_txt,
                      '.mzxml': read_spectrum_from_mzxml}[file_ext]
        except KeyError:
            raise KeyError("Unsupported extension")

        spectrum = reader(spectrum,
                          spectral_intensity_cut_off,
                          precision_digits)
    else:
        spectrum = threshold_n_round_n_aggregate(
            spectrum,
            spectral_intensity_cut_off,
            precision_digits)

    return spectrum
    # spectra = {}
    # spectra['original'] = spectrum
    # spectra['original total intensity'] = sum(spectrum[1])
    #
    # if opt_P:  # cut_off == None
    #     cut_off = remove_lower_quantile(*spectrum, retained_percentage=opt_P)
    #
    # spectra['trimmed'], spectra['left out'] = trim_spectrum(*spectrum,
    #                                                         cut_off=cut_off)
    # spectra['total intensity after trim'] = spectra['trimmed'][1].sum()
    # spectra['trimmed intensity'] = \
    #     spectra['original total intensity'] - \
    #     spectra['total intensity after trim']
    # spectra['cut_off'] = cut_off
    #
    # if _verbose:
    #     print()
    #     print('original total intensity',
    #           spectra['original total intensity'])
    #     print('total intensity after trim',
    #           spectra['total intensity after trim'])
    #     print('trimmed intensity',
    #           spectra['trimmed intensity'])
    #     print()
    # return spectra
