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

from math import fsum
from collections import defaultdict
from six.moves import range, zip
from operator import itemgetter


def cdata2numpyarray(x):
    """Turn c-data into a numpy array.
    Parameters
    ----------
    x : cdata table
        A table of cdata from cffi.

    Returns
    -------
    res : array
        A numpy array of numbers.
    """
    res = np.empty(len(x))
    for i in range(len(x)):
        res[i] = x[i]
    return res


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


def merge_runs(spec1, spec2):
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

# WHAT THE HELL IS THIS?
def aggregate_envelopes(masses, probs, digits=2):
    """Aggregate theoretical envelopes.

    Parameters
    ----------
    masses : array
        An array of isotopologues' masses.
    probs : array
        An array of isotopologues' probabilities.
    digits : int
        The number of significant digits used
        while rounding the masses of isotopologues.

    Returns
    ----------
    out : tuple
        A theoretical spectrum of a given resolution.
    """

    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array([k for k in lists])
    newProbs = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs, op_flags=['readwrite']),
                          newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs


def trim_spectrum(mz, intensity, cut_off=100):
    """Remove peaks below a given cut off.

    Parameters
    ----------
    mz : array
        The m/z values.
    intensity : array
        The intensities to bin.
    cut_off : float
        The cut off value for intensity.

    Returns
    -------
    out : tuple
        A tuple containing the binned experimental spectrum:
        mass over charge values and intensities, both numpy arrays.
    """
    return ((mz[intensity >= cut_off], intensity[intensity >= cut_off]),
            (mz[intensity < cut_off], intensity[intensity < cut_off]))


def round_spectrum(mz, intensity, prec_digits=2):
    """Bin the spectrum by rounding m/z with the same number of significant digits.

    Parameters
    ----------
    mz : array
        The m/z values.
    intensity : array
        The intensities to bin.
    prec_digits : float
        The number of digits after which the floats get rounded.

    Returns
    -------
    out : tuple
        A tuple containing the binned experimental spectrum:
        mass over charge values and intensities, both numpy arrays.
    """
    mz = np.round(mz, prec_digits)
    mz, intensity = aggregate(mz, intensity)
    return mz, intensity


def remove_lower_quantile(mz, intensities, retained_percentage=.95):
    """Remove a portion of the smallest peaks that cover
    1-retained_percentage of the total ion current.

    Parameters
    ----------
    mz : array
        The m/z values.
    intensity : array
        The intensities to bin.
    retained_percentage : float
        The percentage of the original total ion current
        to be retained after trimming.

    Returns
    -------
    effective_cut_off : float
        The minimal peak height.
    """
    mz_res, intensities_res = tuple(
        np.array(x) for x in
        zip(*sorted(zip(mz, intensities), key=itemgetter(1))))

    i = 0
    S = 0.0
    total = intensities_res.sum()
    while True:
        S += intensities_res[i]/total
        if S < 1.0-retained_percentage:
            i += 1
        else:
            break
    effective_cut_off = intensities_res[i]
    return effective_cut_off
