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
from MassTodonPy.Parsers.formula_parser import atomCnt2string


Molecule = namedtuple("Molecule",
                      ("name", "source", "formula", "q", "g"))


class MoleculeMaker(object):
    """A class for obtaining chemical molecules
    that appear in the reaction graph.

    Parameters
    ----------
    precursor : Precursor
        An instance of the precursor class.
    """

    def __init__(self, precursor, blockedFragments=set()):
        self.precursor = precursor
        print(precursor.modifications)
        self.block_fragments(blockedFragments)
        self.make_superatoms()

    def block_fragments(self, blockedFragments):
        self.blockedFragments = set('c0') | blockedFragments
        for i, f in enumerate(self.precursor.fasta):
            if f == 'P':
                self.blockedFragments.add('c' + str(i))
                z_frag_No = len(self.precursor.fasta) - i
                self.blockedFragments.add('z' + str(z_frag_No))

    def __get_AA(self, amino_acid_group, amino_acid_tag, amino_acid_No):
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

    def get_amino_acids(self, direction):
        if direction == "N -> C":
            print('N -> C')
            sA = lCnt({'H': 1})
            for AA_No, AA_tag in enumerate(self.precursor.fasta):
                sA += self.__get_AA('N', AA_tag, AA_No)
                yield sA
                sA = self.__get_AA('C_alpha', AA_tag, AA_No) +\
                    self.__get_AA('C_carbo', AA_tag, AA_No)
            yield sA + lCnt({'O': 1, 'H': 1})
        else:
            N = len(self.precursor.fasta)
            print('C -> N')
            sA = lCnt({'O': 1, 'H': 1})
            for AA_No, AA_tag in enumerate(reversed(self.precursor.fasta)):
                sA += self.__get_AA('C_carbo', AA_tag, AA_No) +\
                      self.__get_AA('C_alpha', AA_tag, N-AA_No)
                yield sA
                sA = self.__get_AA('N', AA_tag, AA_No)
            yield sA + lCnt({'H': 1})

    def make_superatoms(self):
        self.superAtoms = []
        sA = lCnt()
        for AA_No, AA_tag in enumerate(self.precursor.fasta):
            sA += self.__get_AA('N', AA_tag, AA_No)
            self.superAtoms.append(sA)
            sA = self.__get_AA('C_alpha', AA_tag, AA_No) + \
                 self.__get_AA('C_carbo', AA_tag, AA_No)
        sA += lCnt({'O': 1, 'H': 1})
        self.superAtoms.append(sA)
        self.superAtoms[0] += lCnt({'H': 1})

    def protonate(self, frag):
        a, b, c = {'p': (1, 0, 1),
                   'c': (0, -1, 0),
                   'z': (0, 0, 1)}[frag]
        for q in range(1, self.precursor.q + a):
            for g in range(b, self.precursor.q - q + c):
                yield (q, g)

    def unprotonated_molecules(self):
        N = len(self.superAtoms)

        yield ('precursor',
               atomCnt2string(sum(self.superAtoms)),
               len(self.precursor.fasta))

        # Adding one extra hydrogen to meet the definition of a c fragment.
        if 'c' in self.precursor.fragmentation_type:
            cFrag = lCnt({'H': 1})
            for i in range(N-1):
                cFrag += self.superAtoms[i]
                cFrag_tmp = lCnt(cFrag)
                frag_type = 'c' + str(i)
                if frag_type not in self.blockedFragments and not i == 0:
                    yield (frag_type, atomCnt2string(cFrag_tmp), i)

        if 'z' in self.precursor.fragmentation_type:
            zFrag = lCnt()
            for i in range(1, N):
                zFrag += self.superAtoms[N-i]
                zFrag_tmp = lCnt(zFrag)
                frag_type = 'z'+str(i)
                if frag_type not in self.blockedFragments:
                    yield (frag_type, atomCnt2string(zFrag_tmp), i)

    def __iter__(self):
        for molType, formula, sideChainsNo in self.unprotonated_molecules():
            p_c_z = molType[0]
            for q, g in self.protonate(p_c_z):
                potentialChargesNo = \
                    sideChainsNo // self.precursor.distance_charges

                if sideChainsNo % self.precursor.distance_charges > 0:
                    potentialChargesNo += 1
                    # +0000 +0000 00+  at most 3 charges

                if potentialChargesNo >= q:
                    yield Molecule(molType,
                                   self.precursor.name,
                                   formula, q, g)


def make_molecules(precursors):
    """Generate molecules from the Roepstorff Scheme.
    Parameters
    ----------
    precursors
        An iterable full of Precursor objects.

    Returns
    -------
    out : iterator
        An iterator over Molecules generated by the list of precursors.
    """
    for prec in precursors:
        for mol in MoleculeMaker(prec):
            yield mol
