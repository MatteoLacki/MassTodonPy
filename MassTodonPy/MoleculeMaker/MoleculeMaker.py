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
from collections import namedtuple
from MassTodonPy.Data.get_data import get_amino_acids
from MassTodonPy.Parsers.formula_parser import atomCnt2string


class NegativeAtomCount(Exception):
    pass


Molecule = namedtuple("Molecule", ("name", "source", "formula", "q", "g"))


class MoleculeMaker(object):
    """A class for obtaining chemical molecules
    that appear in the reaction graph.

    Parameters
    ----------
    precusor : Precursor
        An instance of the precursor class.
    """

    def __init__(self, precusor):
        self.precusor = precusor
        self.amino_acids = get_amino_acids()
        self.block_fragments()
        self.make_superatoms()
        self.make_unprotonated_molecules()

    def block_fragments(self):
        self.blockedFragments = set()

    def get_amino_acid(self, amino_acid_group, amino_acid_tag, amino_acid_No):
        try:
            amino_acid = \
                self.amino_acids[amino_acid_tag][amino_acid_group] +\
                self.precusor.modifications[amino_acid_No][amino_acid_group]
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

    def make_superatoms(self):
        raise NotImplementedError("Need to know fragmentation type.")

    def protonate(self, frag):
        a, b, c = {'p': (1, 0, 1),
                   'c': (0, -1, 0),
                   'z': (0, 0, 1)}[frag]
        for q in range(1, self.precusor.q + a):
            for g in range(b, self.precusor.q - q + c):
                yield (q, g)

    def make_unprotonated_molecules(self):
        N = len(self.superAtoms)

        # I don't use Precursor here, because this is my kingdom object.
        precursor = ('precursor',
                     atomCnt2string(sum(self.superAtoms)),
                     len(self.precursor.fasta))
        self.unprotonated_molecules = [precursor]

        # c fragments
        # Adding one extra hydrogen to meet the definition of a c fragment.
        cFrag = lCnt({'H': 1})
        for i in range(N-1):
            cFrag += self.superAtoms[i]
            cFrag_tmp = lCnt(cFrag)
            frag_type = 'c' + str(i)
            if frag_type not in self.blockedFragments and not i == 0:
                c_fragment = (frag_type, atomCnt2string(cFrag_tmp), i)
                self.unprotonated_molecules.append(c_fragment)

        # z fragments
        zFrag = lCnt()
        for i in range(1, N):
            zFrag += self.superAtoms[N-i]
            zFrag_tmp = lCnt(zFrag)
            frag_type = 'z'+str(i)
            if frag_type not in self.blockedFragments:
                z_fragment = (frag_type, atomCnt2string(zFrag_tmp), i)
                self.unprotonated_molecules.append(z_fragment)

    def __iter__(self):
        self.molecules = []

        for molType, atomCnt_str, sideChainsNo in self.unprotonated_molecules:
            for q, g in self.protonate(molType[0]):
                potentialChargesNo = \
                    sideChainsNo // self.precusor.distance_charges
                if sideChainsNo % self.precusor.distance_charges > 0:
                    potentialChargesNo += 1
                    # +0000 +0000 00+  at most 3 charges
                if potentialChargesNo >= q:
                    molecule = Molecule(molType, atomCnt_str, q, g)
                    self.molecules.append(molecule)


class CzMoleculeMaker(MoleculeMaker):
    """Enumerate c-z fragments."""

    def block_fragments(self):
        self.blockedFragments = set('c0')
        for i, f in enumerate(self.precusor.fasta):
            if f == 'P':
                self.blockedFragments.add('c' + str(i))
                z_frag_No = len(self.precusor.fasta) - i
                self.blockedFragments.add('z' + str(z_frag_No))

    def make_superatoms(self):
        self.superAtoms = []
        sA = lCnt()
        for amino_acid_No, amino_acid_tag in enumerate(self.precursor.fasta):
            sA += self.get_amino_acid('N',
                                      amino_acid_tag,
                                      amino_acid_No)
            self.superAtoms.append(sA)
            sA = self.get_amino_acid('C_alpha',
                                     amino_acid_tag,
                                     amino_acid_No) + \
                self.get_amino_acid('C_carbo',
                                    amino_acid_tag,
                                    amino_acid_No)
        sA += lCnt({'O': 1, 'H': 1})
        self.superAtoms.append(sA)
        self.superAtoms[0] += lCnt({'H': 1})


# TODO: add other fragmentation schemes
def get_molecules(precursors):
    """Generate molecules from the Roepstorff Scheme.
    Parameters
    ----------
    precursors
        An iterable full of Precursor objects.

    Returns
    -------
    out : class
        A list of molecules.

    Warning
    -------
        At present only the 'cz' fragmentation is supported.
    """

    F = [ mol for mole in CzMoleculeMaker(precusor) for precusor in precursors]

    return F.molecules
