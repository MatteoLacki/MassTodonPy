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

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator


class Molecule(object):
    """A class representing reaction products.

    Parameters
    ==========
    name : string
        The conventional name of the molecule, e.g. precursor, c11, z23, a12.
    source : Precursor
        An instance of the Precursor object: the parent ion for this molecule.
    formula : Formula
        An instance of the Formula, a chemical formula object.
    q : int
        Charge of the reaction product.
    g : int
        Quenched charge of the reaction product.
    """
    iso_calc = IsotopeCalculator(joint_probability=.999,
                                 mz_digits=infinity)

    @classmethod
    def reset_isotope_calculator(cls, **args):
        """Reset the IsotopeCalculator with its usual initial arguments."""
        cls.iso_calc = IsotopeCalculator(**args)

    def __init__(self, name, source, formula, q=0, g=0):
        self.name = name
        self.source = source
        self.formula = formula
        self.q = int(q)
        self.g = int(g)
        self.intensity = 0.0

    #TODO generalize to abxy
    def _molType_position_cleavageSite(self):
        mt = self.name[0]
        if mt == 'p':
            return None
        else:
            po = int(self.name[1:])
            fasta_len = len(self.source.fasta)
            cs = None if mt == 'p' else \
                   po if mt == 'c' else fasta_len - po
            return mt, po, cs

    @property
    def monoisotopic_mz(self):
        return self.iso_calc.get_monoisotopic_mz(self.formula, self.q, self.g)

    @property
    def mean_mz(self):
        return self.iso_calc.get_mean_mz(self.formula, self.q, self.g)

    def isotopologues(self,
                      joint_probability=.999,
                      mz_digits=infinity,
                      **kwds):
        return self.iso_calc.get_envelope(self.formula,
                                          joint_probability,
                                          self.q,
                                          self.g,
                                          mz_digits,
                                          memoize=True)

    def __repr__(self):
        return "({name} {source.name} q={q} g={g} I={I_int})".format(
            I_int=int(self.intensity),
            **self.__dict__)

    def __hash__(self):
        return hash((self.name,
                     self.source,
                     str(self.formula),
                     self.q,
                     self.g))
# Python 3 does not use cmp anymore.
