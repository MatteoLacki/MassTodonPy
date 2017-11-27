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

from six.moves import range
from MassTodonPy.Formula.Formula import Formula
from collections import namedtuple


Molecule = namedtuple("Molecule",
                      ("name", "source", "formula", "q", "g"))


class MoleculeMaker(object):
    """
    A class for obtaining chemical molecules
    that appear in the reaction graph.
    """
    def __init__(self,
                 precursor,
                 blocked_fragments=set(['c0']),
                 fragmentation_type="cz",
                 distance_charges=5):
        self.precursor = precursor
        self.blocked_fragments = blocked_fragments
        self.fragmentation_type = fragmentation_type
        self.distance_charges = distance_charges
        for i, f in enumerate(self.precursor.fasta):
            if f == 'P':
                self.blocked_fragments.add('c' + str(i))
                z_frag_No = len(self.precursor.fasta) - i
                self.blocked_fragments.add('z' + str(z_frag_No))

    def protonate(self, frag):
        a, b, c = {'p': (1, 0, 1),
                   'c': (0, -1, 0),
                   'z': (0, 0, 1)}[frag]
        for q in range(1, self.precursor.q + a):
            for g in range(b, self.precursor.q - q + c):
                yield (q, g)

    def c_fragments(self):
        """Generate c fragments."""
        # extra hydrogen to meet
        # the definition of a c fragment.
        # not as H on the N terminus!
        # H on the N terminus is added
        # in the precursor __getitem__ method
        atom_cnt = Formula('H1')
        for aa_cnt in range(len(self.precursor.fasta)):
            atom_cnt += self.precursor[aa_cnt, 'N']
            name = 'c' + str(aa_cnt)
            if name not in self.blocked_fragments:
                yield (name, str(atom_cnt))
            atom_cnt += self.precursor[aa_cnt, 'C_alpha']
            atom_cnt += self.precursor[aa_cnt, 'C_carbo']

    def z_fragments(self):
        """Generate z fragments."""
        atom_cnt = Formula()
        for aa_cnt in range(len(self.precursor.fasta) - 1, -1, -1):
            atom_cnt += self.precursor[aa_cnt, 'C_carbo']
            atom_cnt += self.precursor[aa_cnt, 'C_alpha']
            side_chain_len = len(self.precursor.fasta) - aa_cnt
            name = 'z' + str(side_chain_len)
            if name not in self.blocked_fragments:
                yield (name, str(atom_cnt))
            atom_cnt += self.precursor[aa_cnt, 'N']

    def uncharged_molecules(self):
        """Generate uncharged molecules."""
        yield ('precursor',  # name
               str(self.precursor.formula))
        if 'c' in self.fragmentation_type:
            for c_frags in self.c_fragments():
                yield c_frags
        if 'z' in self.fragmentation_type:
            for z_frags in self.z_fragments():
                yield z_frags

    def charged_molecules(self):
        """Generate charged molecules."""
        for name, atom_cnt in self.uncharged_molecules():
            p_c_z = name[0]
            if p_c_z is not 'p':
                side_chain_len = int(name[1:])
            else:
                side_chain_len = len(self.precursor.fasta)
            for q, g in self.protonate(p_c_z):
                potential_charges_cnt = \
                    side_chain_len // self.distance_charges
                if side_chain_len % self.distance_charges > 0:
                    potential_charges_cnt += 1
                    # +0000 +0000 00+  at most 3 charges
                if potential_charges_cnt >= q:
                    yield Molecule(name,
                                   self.precursor.name,
                                   atom_cnt, q, g)


def get_molecules(precursors,
                  blocked_fragments=set(['c0']),
                  fragmentation_type="cz",
                  distance_charges=5):
    """
    Generate molecules from the Roepstorff Scheme.

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
        mol_maker = MoleculeMaker(prec,
                                  blocked_fragments,
                                  fragmentation_type,
                                  distance_charges)

        for mol in mol_maker.charged_molecules():
            yield mol
