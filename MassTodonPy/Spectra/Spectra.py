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
from collections import namedtuple, Counter
from MassTodonPy.Spectra.operations import aggregate

TheoreticalSpectrum = namedtuple('TheoreticalSpectrum', 'mz probability')
ExperimentalSpectrum = namedtuple('ExperimentalSpectrum', 'mz intensity')


def threshold_n_round_n_aggregate(spectrum,
                                  spectral_intensity_cut_off=0.0,
                                  precision_digits=2):
    """
    Apply the peak intensity thresholding, rounding m/z values,
    and aggregate the spectrum.
    Parameters
    ----------
    path : str
        Path to the *.txt text file containing the mass spectrum.
    spectral_intensity_cut_off : float
        The cut off value for peak intensity.
    precision_digits : float
        The number of digits after which the floats get rounded.
    Returns
    -------
    out : Spectrum
    """
    mzs, intensities = spectrum
    mzs = mzs[
        intensities >= spectral_intensity_cut_off]

    intensities = intensities[
        intensities >= spectral_intensity_cut_off]

    mzs = np.round(mzs, precision_digits)
    mzs, intensities = aggregate(mzs, intensities)
    return ExperimentalSpectrum(mz=mzs, intensity=intensities)


def threshold_round_and_aggregate_peaks(peaks,
                                        spectral_intensity_cut_off=0.0,
                                        precision_digits=2):
    """
    Apply thresholding on peaks,
    round their m/z values,
    and aggregate the result.

    Parameters
    =======
    peaks : iterable
        Iterates over tuples (m_over_z, intensity).
    spectral_intensity_cut_off : float

    precision_digits : float
        The number of digits after which the floats get rounded.
    Remarks
    =======
    This function is purely pythonic (if the peak iterable is).
    It is 8 times slower than the numpy version.
    """
    nice_spectrum = Counter()
    for mz, intensity in peaks:
        if intensity >= spectral_intensity_cut_off:
            nice_spectrum[round(mz, precision_digits)] += intensity
    return ExperimentalSpectrum(
            mass=np.array(list(nice_spectrum.keys())),
            intensity=np.array(list(nice_spectrum.values())))


#
# def read_spectrum_from_mzxml(path,
#                              spectral_intensity_cut_off=0.0,
#                              precision_digits=2,
#                              spectra_reader=get_data_from_mzxml_faster):
#     """
#     Read spectrum form an mzXml and merge runs of the instrument.
#     Parameters
#     ----------
#     path : str
#         Path to the mzXml file containing the mass spectrum.
#     spectral_intensity_cut_off : float
#         The cut off value for peak intensity.
#     precision_digits : float
#         The number of digits after which the floats get rounded.
#     Returns
#     -------
#     out : Spectrum
#     """
#     assert isinstance(spectral_intensity_cut_off, (int, float)) and\
#         spectral_intensity_cut_off >= 0,\
#         'Intensity cut off should be nonnegative.'
#
#     assert isinstance(precision_digits, int) and precision_digits >= 0,\
#         'Precision digits should be a natural number or zero.'
#
#     mz, intensity = reduce(
#         stack_measures,
#         spectra_reader(path,
#                        spectral_intensity_cut_off,
#                        precision_digits))
#     return Spectrum(mz=mz, intensity=intensity)
