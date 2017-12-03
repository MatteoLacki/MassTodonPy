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
    iso_calc = IsotopeCalculator(mz_precision=3,
                                 joint_probability=.999)

    @classmethod
    def reset_isotope_calculator(cls, **args):
        """Reset the IsotopeCalculator with its usual initial arguments."""
        cls.iso_calc = IsotopeCalculator(**args)

    def __init__(self, name, source, formula, q=0, g=0):
        self.name = name
        self.source = source
        self.formula = formula
        self.q = q
        self.g = g

    @property
    def monoisotopic_mz(self):
        return self.iso_calc.get_monoisotopic_mz(self.formula, self.q, self.g)

    @property
    def mean_mz(self):
        return self.iso_calc.get_mean_mz(self.formula, self.q, self.g)

    def isotopologues(self,
                      joint_probability=.999,
                      mz_precision=3):
        return self.iso_calc.get_envelope(self.formula,
                                          joint_probability,
                                          self.q,
                                          self.g,
                                          mz_precision,
                                          memoize=True)

    def __repr__(self):
        out = "Molecule {name} out of {source}:\n".format(**self.__dict__)
        out += "\t{}\n".format(self.formula.__repr__())
        out += "\tCharge = {q}\n\tQuenched charge = {g}\n".format(**self.__dict__)
        return out

    def __hash__(self):
        return hash((self.name,
                     self.source,
                     str(self.formula),
                     self.q,
                     self.g))
# Python 3 does not use cmp anymore.
