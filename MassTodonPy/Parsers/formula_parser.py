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
from collections import Counter
from MassTodonPy.Data.get_isotopes import get_elements


def get_formula_parser(pattern='([A-Z][a-z]?)([0-9]*)'):
    """Set up the formula.

    Parameters
    ----------
    pattern : str
        A regular expression that describes elements.

    Returns : function
        A formula parser.
    """

    pattern = re.compile(pattern)

    def parse(atomCnt_str):
        """Parses chemical formula based on the class pattern definition.

        Parameters
        ----------
        atomCnt_str : str
            The chemical formula string.

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
        atomCnt = Counter()
        for elemTag, cnt in re.findall(pattern, atomCnt_str):
            if elemTag in get_elements():
                if cnt == '':
                    cnt = 1
                else:
                    cnt = int(cnt)
                atomCnt[elemTag] += cnt
            else:
                raise AttributeError
        return atomCnt

    return parse


parse_formula = get_formula_parser()


def atom_cnt_2_string(atomCnt):
    """Translate a dictionary of atom counts into a uniquely defined string.

    Parameters
    ----------
    atomCnt : Counter
        The chemical formula counter.

    Returns
    -------
    out : str
        A chemical formula string.
    """
    keys = atomCnt.keys()
    return "".join(el+str(atomCnt[el]) for el in sorted(keys))
