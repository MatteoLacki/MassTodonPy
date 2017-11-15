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

from collections import defaultdict
from linearCounter.linearCounter import linearCounter as lCnt
from MassTodonPy.Data.get_amino_acids import get_amino_acids


class NegativeAtomCount(Exception):
    pass


class Precursor(object):
    """A class for storing information on precursors.

    fasta : str
        The fasta of the studied molecular species.

    Q : int
        The charge of the precursor ion.

    distance_charges : int
        The minimal distance between charges on the protein.
        If set to 5, at least 4 amino acids must lay between
        consecutive *charged* amino acids.

    modifications : dictionary
        A dictionary of modifications of amino acids.

        Key: amino acid number in fasta sequence
             (beginning with N terminus, fininshing with C terminus).

        Value: a dictionary with group modifications.
            Keys : C_carbo, C_alpha, or N.
            Value: atom count in form of a linearCounter.
    """
    __slots__ = ("name", "fasta", "q",
                 "fragmentation_type",
                 "distance_charges",
                 "modifications",
                 "modifications2")

    amino_acids = get_amino_acids()

    def __init__(self, name, fasta, q,
                 fragmentation_type="cz",
                 distance_charges=5,
                 modifications={}):
        self.name = name
        self.fasta = fasta
        self.q = q
        self.fragmentation_type = fragmentation_type
        self.distance_charges = distance_charges

        # turning modifications into linear counters
        self.modifications2 = defaultdict(
            lambda: defaultdict(lambda: lCnt()),
            {k - 1: defaultdict(lambda: lCnt(),
                                {name: lCnt(atomCnt)
                                 for name, atomCnt in v.items()})
             for k, v in modifications.items()})

        self.modifications = {(k - 1, name): lCnt(atomCnt)
                              for k, v in modifications.items()
                              for name, atomCnt in v.items()}

    def __getitem__(self, amino_acid_group, amino_acid_tag, amino_acid_No):
        try:
            amino_acid = \
                self.amino_acids[amino_acid_tag][amino_acid_group] +\
                self.precursor.modifications[amino_acid_No][amino_acid_group]
            print(self.precursor.modifications)

        except KeyError:
            print(amino_acid_No, amino_acid_tag, amino_acid_group)
            print(self.amino_acids[amino_acid_tag][amino_acid_group])
        try:
            if any(count < 0 for element, count in amino_acid.items()):
                raise NegativeAtomCount("Attention: your modification had an unexpected effect.\
                Part of your molecule now has negative atom count.\
                Bear that in mind while publishing your results.")
        except UnboundLocalError:
            amino_acid = lCnt()
        return amino_acid
