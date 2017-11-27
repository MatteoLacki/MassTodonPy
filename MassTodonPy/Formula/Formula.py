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

from __future__ import absolute_import, division, print_function
from future.builtins import super

from MassTodonPy.Formula.Parse import get_formula_parser
from MassTodonPy.Formula.LinearDict import LinearDict


class Formula(LinearDict):
    parse = get_formula_parser()

    def reset_parser(self, pattern='([A-Z][a-z]?)([0-9]*)'):
        self.__class__.parse = get_formula_parser(pattern)

    def __init__(self, formula={}):
        if isinstance(formula, str):
            formula = self.__class__.parse(formula)
        super().__init__(formula)

    def __str__(self):
        return "".join(element + str(count)
                       for element, count in
                       sorted(self._storage.items()))

    def __repr__(self):
        out = self._storage.__repr__()
        out = "Formula({})".format(out[1:-1])
        return out
