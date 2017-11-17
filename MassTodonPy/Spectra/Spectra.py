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
from six.moves import reduce
from pyteomics import mzxml  # >= 3.41

Spectrum = namedtuple('Spectrum', 'mz intensity')
TheoreticalSpectrum = namedtuple('TheoreticalSpectrum', 'mz probability')


def aggregate(keys, values=None):
    """Aggregate values with the same keys.

    Parameters
    ----------
    keys : array
        Keys, usually m/z values.
    values : array
        Values to aggregate, usually intensities.

    Returns
    -------
    out : tuple
        A tuple containing unique keys and aggregated values.
    """
    uniqueKeys, indices = np.unique(keys, return_inverse=True)
    return uniqueKeys, np.bincount(indices, weights=values)


def stack_spectra(spec1, spec2):
    """Merge two spectra into one, aggregating out the common m/z values.

    Parameters
    ----------
    spec1 : tuple of two numpy arrays
        A mass spectrum:
        an array of m/z ratios and an array of corresponding intensities.
    spec2 : tuple of two numpy arrays
        A mass spectrum:
        an array of m/z ratios and an array of corresponding intensities.
    """
    masses_over_charge = np.concatenate((spec1[0], spec2[0]))
    intensities = np.concatenate((spec1[1], spec2[1]))
    return aggregate(masses_over_charge, intensities)


def get_distributions_from_mzxml(path,
                                 spectral_intensity_cut_off=0.0,
                                 precision_digits=2):
    """
    Generate a sequence of rounded and trimmed spectra from
    individual runs of the instrument.

    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    spectral_intensity_cut_off : float
        The cut off value for peak intensity.
    precision_digits : float
        The number of digits after which the floats get rounded.
    Returns
    -------
    out : generator
        Generates tuples of numpy arrays corresponding to different runs
        of the experimental spectrum.
    """
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mzs = spectrum['m/z array']

            intensities = spectrum['intensity array']

            mzs = mzs[intensities >=
                      spectral_intensity_cut_off]

            mzs = np.round(mzs, precision_digits)
            intensities = intensities[intensities >=
                                      spectral_intensity_cut_off]
            yield mzs, intensities


def read_spectrum_from_mzxml(path,
                             spectral_intensity_cut_off=0.0,
                             precision_digits=2):
    """
    Read spectrum form an mzXml and merge runs of the instrument.
    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    spectral_intensity_cut_off : float
        The cut off value for peak intensity.
    precision_digits : float
        The number of digits after which the floats get rounded.
    Returns
    -------
    out : Spectrum
    """
    mz, intensity = reduce(
        stack_spectra,
        get_distributions_from_mzxml(path,
                                     spectral_intensity_cut_off,
                                     precision_digits))
    return Spectrum(mz=mz, intensity=intensity)


def read_spectrum_from_txt(path,
                           spectral_intensity_cut_off=0.0,
                           precision_digits=2):
    """
    Read spectrum from a text file.
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
    mzs = []
    intensities = []
    with open(path) as f:
        for line in f:
            line = line.split()
            intensity = float(line[1])
            if intensity >= spectral_intensity_cut_off:
                mzs.append(float(line[0]))
                intensities.append(intensity)
    mzs = np.round(np.array(mzs), precision_digits)
    intensities = np.array(intensities)
    mzs, intensities = aggregate(mzs, intensities)
    return Spectrum(mz=mzs, intensity=intensities)


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

    return Spectrum(mz=mzs, intensity=intensities)


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
    return Spectrum(mass=np.array(list(nice_spectrum.keys())),
                    intensity=np.array(list(nice_spectrum.values())))
