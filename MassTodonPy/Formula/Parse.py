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
import re
from MassTodonPy.Data.get_isotopes import get_elements


def get_pattern(pattern='([A-Z][a-z]?)([0-9]*)'):
    return re.compile(pattern)


def parse(formula, pattern):
    """Parses chemical formula based on the class pattern definition.

    Parameters
    ----------
    formula : str
        The chemical formula string.
    pattern : str
        The 'compiled' pattern for parsing chemical formulas.

    Returns
    -------
    atomCnt : Counter
        A counter with elements for keys and atom counts for values.

    Examples
    --------
        >>> FP = formulaParser()
        >>> FP('C100H202')
        Counter({'C': 100, 'H': 202})
    """
    atomCnt = {}
    for elemTag, cnt in re.findall(pattern, formula):
        if elemTag in get_elements():
            if cnt == '':
                cnt = 1
            else:
                cnt = int(cnt)
            atomCnt[elemTag] = cnt
        else:
            raise AttributeError
    return atomCnt
