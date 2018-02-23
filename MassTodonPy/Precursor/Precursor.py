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

from MassTodonPy.Data.get_amino_acids import get_amino_acids
from MassTodonPy.Formula.Formula import Formula
from MassTodonPy.Molecule.Molecule import Molecule


def flatten_modification(mod):
    return ''.join(''.join((str(n), t, str(v)))
                   for (n, t), v in mod.items())


class Precursor(object):
    """Make precursor.

    Parameters
    ==========
    fasta : str
        The fasta of the studied molecular species.
    name : str
        The name of the precursor molecule.
    charge : int
        The charge of the precursor ion.
    modifications : dictionary
        A dictionary of modifications of amino acids.
        Key: amino acid number in fasta sequence
             (from N to C termini).
        Value: a dictionary with group modifications.
            Keys : C_carbo, C_alpha, or N.
            Value: atom count in form of a linearCounter.
    fragments: str
        For now 'cz' only, but we are working on it.
    blocked_fragments : list
        Fragments you don't want to include, e.g. 'z5'.
    block_prolines : boolean
        Should we block prolines?
    distance_charges :
        The minimal distance between charges on the fasta sequence.
        Defaults to charges being 4 amino acids apart.
    kwds :
        Settings for other methods.
    """
    amino_acids = get_amino_acids()  # residues only!

    def __init__(self,
                 fasta,
                 charge,
                 name="",
                 modifications={},
                 fragments="cz",
                 blocked_fragments=set(['c0']),
                 block_prolines=True,
                 distance_charges=5,
                 **kwds):
        self.name = name
        self.fasta = fasta
        self.q = int(charge)
        self.groups = ('N', 'C_alpha', 'C_carbo')
        self.group2frag = dict(start=dict(N='y', C_alpha='z', C_carbo='x'),
                               end=dict(N='c', C_alpha='a', C_carbo='b'))
        self.modifications = {(number - 1, group): Formula(atom_cnt)
                              for number, mods in modifications.items()
                              for group, atom_cnt in mods.items()}
        self.formula = sum(self[number, group]
                           for number in range(len(self))
                           for group in self.groups)
        self.blocked_fragments = blocked_fragments
        self.fragments = fragments
        self.distance_charges = int(distance_charges)
        if block_prolines:
            for i, f in enumerate(self.fasta):
                if f == 'P':
                    self.blocked_fragments.add('c' + str(i))
                    z_frag_No = len(self) - i
                    self.blocked_fragments.add('z' + str(z_frag_No))

    def _get_amino_acid(self, number, group):
        """Get amino acid of the precursor."""
        formula = self.amino_acids[(self.fasta[number], group)] +\
                  self.modifications.get((number, group), 0)
        # Modifying termini: Kaltashov O., Mass Spectrometry in Biophysics
        if number is 0 and group is 'N':
            formula['H'] += 1  #  H for N terminus
        if number is len(self)-1 and group is 'C_carbo':
            formula['O'] += 1  # additional H and O
            formula['H'] += 1  # for the C terminus
        formula.check_positivity()
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
        out = "({name} {fasta} q={q}".format(**self.__dict__)
        out += ' modified)' if self.modifications else ')'
        return out

    def __len__(self):
        return len(self.fasta)

    def _protonate(self, frag):
        a, b, c = {'p': (1, 0, 1),
                   'c': (0, -1, 0),
                   'z': (0, 0, 1)}[frag]
        for q in range(1, self.q + a):
            for g in range(b, self.q - q + c):
                yield (q, g)

    def a_fragments(self):
        """Generate a fragments."""
        pass

    def b_fragments(self):
        """Generate b fragments."""
        pass

    def c_fragments(self):
        """Generate c fragments."""
        # 'H1' to be a 'c' fragment, not the H on the N terminus
        formula = Formula('H1')
        for number in range(len(self)):
            formula += self[number, 'N']
            name = 'c' + str(number)
            if name not in self.blocked_fragments:
                yield (name, formula.copy())
            formula += self[number, 'C_alpha']
            formula += self[number, 'C_carbo']

    def x_fragments(self):
        """Generate x fragments."""
        pass

    def y_fragments(self):
        """Generate y fragments."""
        pass

    def z_fragments(self):
        """Generate z fragments."""
        formula = Formula()
        for number in range(len(self)-1, -1, -1):
            formula += self[number, 'C_carbo']
            formula += self[number, 'C_alpha']
            side_chain_len = len(self) - number
            name = 'z' + str(side_chain_len)
            if name not in self.blocked_fragments:
                yield (name, formula.copy())
            formula += self[number, 'N']

    def uncharged_molecules(self):
        """Generate uncharged molecules."""
        yield ('precursor', self.formula.copy())
        for frag_type in self.fragments:
            frags = getattr(self, frag_type + '_fragments')()
            for frag in frags:
                yield frag

    def molecules(self):
        """Generate charged molecules."""
        for name, formula in self.uncharged_molecules():
            if name[0] is not 'p':
                side_chain_len = int(name[1:])
            else:
                side_chain_len = len(self)
            #TODO: derive rules for protonation of fragments other than cz
            for q, g in self._protonate(name[0]):
                potential_charges_cnt = \
                    side_chain_len // self.distance_charges
                if side_chain_len % self.distance_charges > 0:
                    potential_charges_cnt += 1
                    # +0000 +0000 00+  at most 3 charges
                if potential_charges_cnt >= q:
                    yield Molecule(name, self, formula, q, g)

    def __hash__(self):
        """Get a hash from the precursor's unique id.

        Unique id consists of a name, fasta, charge, and modifications.
        """
        return hash((self.name,
                     self.fasta,
                     self.q, flatten_modification(self.modifications)))
