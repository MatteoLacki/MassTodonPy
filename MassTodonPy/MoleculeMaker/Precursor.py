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
    __slots__ = ("name", "fasta", "q", "modifications", "atom_cnt")

    amino_acids = get_amino_acids()
    empty_count = lCnt()

    def __init__(self, name, fasta, q, modifications={}):
        self.name = name
        self.fasta = fasta
        self.q = q
        self.modifications = {
            (aa_no - 1, group_name): lCnt(atom_cnt)
            for aa_no, mods in modifications.items()
            for group_name, atom_cnt in mods.items()}
        self.atom_cnt = sum(self.get_AA(i, g)
                            for i in range(len(fasta))
                            for g in ('N', 'C_alpha', 'C_carbo'))

    def get_AA(self, amino_acid_No, amino_acid_group):
        assert amino_acid_group in ('N', 'C_alpha', 'C_carbo'),\
            "The required group is not in ('N', 'C_alpha', 'C_carbo')"

        assert isinstance(amino_acid_No, int) and \
            0 <= amino_acid_No <= len(self.fasta),\
            "The required amino acid number is invalid."

        amino_acid = self.fasta[amino_acid_No]
        atom_cnt = \
            self.amino_acids[(amino_acid, amino_acid_group)] +\
            self.modifications.get((amino_acid_No,
                                    amino_acid_group),
                                   self.empty_count)

        # self.amino_acids contains only amino acid's residues.
        # these need to be modified on the C and N termini
        # Ref.: Kaltashov O. Eyles S.J., Mass Spectrometry in Biophysics
        if amino_acid_No == 0 and amino_acid_group == 'N':
            atom_cnt['H'] += 1  # the additional H for N terminus
        if amino_acid_No == len(self.fasta) - 1 and\
           amino_acid_group == 'C_carbo':
            atom_cnt['O'] += 1  # the additional 0H
            atom_cnt['H'] += 1  # for the C terminus

        if any(count < 0 for element, count in atom_cnt.items()):
            raise NegativeAtomCount("Attention: your modification had an unexpected effect.\
            Part of your molecule now has negative atom count.\
            Good wishes on trying to publish your results.")

        return atom_cnt

    def __getitem__(self, aa_no):
        return sum(self.get_AA(aa_no, aa_g)
                   for aa_g in ('N', 'C_alpha', 'C_carbo'))
