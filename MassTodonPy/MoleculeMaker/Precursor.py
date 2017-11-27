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

from MassTodonPy.Data.get_amino_acids import get_amino_acids
from MassTodonPy.Formula.Formula import Formula


class NegativeAtomCount(Exception):
    pass


class Precursor(object):
    """A class for storing information on precursors.
    name : str
        The name of the precursor molecule.
    fasta : str
        The fasta of the studied molecular species.
    q : int
        The charge of the precursor ion.
    modifications : dictionary
        A dictionary of modifications of amino acids.

        Key: amino acid number in fasta sequence
             (beginning with N terminus, fininshing with C terminus).

        Value: a dictionary with group modifications.
            Keys : C_carbo, C_alpha, or N.
            Value: atom count in form of a linearCounter.
    """
    amino_acids = get_amino_acids()

    def __init__(self, name, fasta, q, modifications={}):
        self.name = name
        self.fasta = fasta
        self.q = q
        self.groups = ('N', 'C_alpha', 'C_carbo')
        self.modifications = {(number - 1, group): Formula(atom_cnt)
                              for number, mods in modifications.items()
                              for group, atom_cnt in mods.items()}
        self.formula = sum(self[number, group]
                           for number in range(len(fasta))
                           for group in self.groups)

    def _get_amino_acid(self, number, group):
        """Get amino acid of the precursor."""
        amino_acid = self.fasta[number]
        formula = self.amino_acids[(amino_acid, group)] +\
                  self.modifications.get((number, group), 0)
        # self.amino_acids contains only amino acid's residues.
        # these need to be modified on the C and N termini
        # Ref.: Kaltashov O. Eyles S.J., Mass Spectrometry in Biophysics
        if number is 0 and group is 'N':
            formula['H'] += 1  # the additional H for N terminus
        if number is len(self.fasta)-1 and group is 'C_carbo':
            formula['O'] += 1  # the additional 0H
            formula['H'] += 1  # for the C terminus
        if any(count < 0 for element, count in formula.items()):
            raise NegativeAtomCount("Attention: your modification had an unexpected effect.\
            Part of your molecule now has negative atom count.\
            Good wishes on trying to publish your results.")
        return formula

    def __getitem__(self, key):
        """Get amino acid of the precursor.

        Parameters
        ==========
        key : int or tuple(int, str)
            The string should be C_carbo, C_alpha, or N.
            The int describes the number of the amino acid, starting from
            zero, counting from N terminus to C terminus.
        Returns
        =======
        out : Formula
        """
        try:
            return self._get_amino_acid(*key)
        except TypeError:
            return sum(self._get_amino_acid(key, group)
                       for group in self.groups)
        except ValueError:
            raise KeyError("Supply '(number, group)' or just 'group'.")

    def __repr__(self):
        out = "Precursor:\n\tname:\t{}\n".format(self.name)
        out += "\tfasta:\t{}\n".format(self.fasta)
        out += "\t{}\n".format(self.formula.__repr__())
        out += "\tcharge:\t{}\n".format(self.q)
        out += "\tmodifications:\t{}".format(self.modifications.__repr__())
        return out
