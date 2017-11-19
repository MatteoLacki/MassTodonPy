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
from six.moves import zip
from operator import itemgetter


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
