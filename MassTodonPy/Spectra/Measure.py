from bisect import bisect_left
import numpy as np
from operator import itemgetter
from six.moves import zip

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Misc.strings import repr_long_list

class Measure(object):
    """Store a discrete finite measure with atoms in R."""

    def __init__(self, atoms=np.array([]), masses=np.array([])):
        """Initialize a measure.

        Parameters
        ----------
        atoms : numpy array
            The atoms upon which the measure holds the mass.
        masses : numpy array
            The masses on atoms.

        """
        self.atoms = np.array(atoms)
        self.masses = np.array(masses)

    def __has_type_of(self, other):
        """Assert that 'self' and 'other' have the same type."""
        assert self.__class__.__name__ is other.__class__.__name__,\
             "\tIllegal to add class {0} to class {1}.\n".format(
                other.__class__.__name__,\
                self.__class__.__name__)

    def __add__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        atoms = np.concatenate((self.atoms, other.atoms))
        masses = np.concatenate((self.masses, other.masses))
        self.__has_type_of(other)
        new_measure = self.__class__(atoms, masses)
        new_measure.__aggregate()
        return new_measure

    def __radd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        if other is 0:
            return self
        else:
            return self.__add__(other)

    def __iadd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        self.__has_type_of(other)
        self.atoms = np.concatenate((self.atoms, other.atoms))
        self.masses = np.concatenate((self.masses, other.masses))
        self.__aggregate()
        return self

    def __aggregate(self):
        """Aggregate masses with the same atoms."""
        self.atoms, indices = np.unique(self.atoms, return_inverse=True)
        self.masses = np.bincount(indices, weights=self.masses)

    def round_atoms(self, precision=infinity):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
            Defaults to 'inf', which prevents any rounding.

        """
        if precision is not infinity:
            self.atoms = np.around(self.atoms, precision)
            self.__aggregate()

    def trim(self, cut_off):
        """Trim masses below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        if cut_off > 0:
            self.atoms = self.atoms[self.masses >= cut_off]
            self.masses = self.masses[self.masses >= cut_off]

    def split_measure(self, cut_off):
        """Split measure into two according to the cut off on masses.

        Retain the measure with masses greater or equal to the cut off.
        Parameters
        ----------
        cut_off : float
        Returns
        ----------
        other : Measure
            A measure with masses strictly below the cut off.

        """
        other = self.__class__(self.atoms[self.masses < cut_off],
                               self.masses[self.masses < cut_off])
        self.trim(cut_off)
        return other

    def get_P_set_cut_off(self, P=.99):
        """Get the cut off resulting in optimal P-set.

        Parameters
        ----------
        P : float
            The percentage of the initial total value that the new measure
            will contain. The new measure contains only atoms with
            heighest masses.

        """
        assert 0.0 <= P <= 1.0, "Wrong P for P-optimal set."
        total_value = self.masses.sum()
        i = 0
        S = 0.0
        for intensity in sorted(self.masses):
            S += intensity/total_value
            if S < P:
                break
        return intensity

    def __repr__(self):
        """Represent the measure."""
        return "{0}:\n\t* {1}\n\t* {2}\n".format(self.__class__.__name__,
                                                 repr_long_list(self.atoms),
                                                 repr_long_list(self.masses))

    def __len__(self):
        """Get size of the measure: the number of atoms."""
        return len(self.atoms)

    def __iter__(self):
        """Iterate over pairs (atom, mass)."""
        return zip(self.atoms, self.masses)

    def __getitem__(self, L_R):
        """Filter atoms between 'L' and 'R', where 'L_R=(L, R)'."""
        try:
            L, R = L_R
            idx = bisect_left(self.atoms, L)
            while self.atoms[idx] <= R:
                yield self.atoms[idx], self.masses[idx]
                idx += 1
        except IndexError:
            return

    def total_mass(self):
        return self.masses.sum()
