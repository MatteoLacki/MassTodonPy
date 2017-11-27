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
from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator

class Molecule(object):
    isotope_calculator = IsotopeCalculator()

    def __init__(self, name, source, formula, q=0, g=0):
        self.name = name
        self.source = source
        self.formula = formula
        self.q = q
        self.g = g

    def reset_isotope_calculator(self):
        pass

    @property
    def monoisotopic_mz(self):
        pass

    @property
    def mean_mz(self):
        pass

    def isotopologues(self):
        pass

    def __repr__(self):
        out = "Molecule {name} out of {source}:\n".format(**self.__dict__)
        out += "\t{}\n".format(self.formula.__repr__())
        out += "\tCharge = {q}\n\tQuenched charge = {g}\n".format(**self.__dict__)
        return out
