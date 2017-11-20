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
from six.moves import zip
from operator import itemgetter


def parse_path(path):
    """
    Parse path to the file.
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


def aggregate(keys, values=None):
    """
    Aggregate values with the same keys.
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


def stack_measures(measure1, measure2):
    """
    Stack two discrete atomic measures on top of each other.
    Parameters
    ----------
    measure1, measure2 : tuples
        A measure is a tuple.
        The first coordinate contains a numpy array with the support.
        The second - an array with the measure values.
    """
    new_support = np.concatenate((measure1[0], measure2[0]))
    new_values = np.concatenate((measure1[1], measure2[1]))
    return aggregate(new_support, new_values)


def round_n_trim(support,
                 values,
                 support_precision=2,
                 value_cut_off=0.0):
    """
    Round the measure's support to a given precision
    and trim the values lower than cut off value.
    Parameters
    ----------
    support : numpy array
        The atoms of the measure.
    value : numpy array
        The values of the measure.
    value_cut_off : float
        The cut off for measure's values.
    support_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    Returns
    -------
    out : a measure tuple
    """
    support = support[values >= value_cut_off]
    support = np.round(support, support_precision)
    values = values[values >= value_cut_off]
    return support, values


def retained_intensity(mzs, intensities, retained_percentage=.95):
    """
    Get the intensity threshold that retains a given percentage of
    the joint intensity of the given spectrum.

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
        zip(*sorted(zip(mzs, intensities), key=itemgetter(1))))

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