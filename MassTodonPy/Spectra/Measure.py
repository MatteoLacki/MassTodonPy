import numpy as np
from operator import itemgetter
from six.moves import zip

class Measure(object):
    """Store a discrete finite measure with support in R."""

    def __init__(self, support=np.array(), values=np.array()):
        """Initialize a measure.

        Parameters
        ----------
        support : a numpy array
        """
        self.__support = np.array(support)
        self.__values = np.array(values)

    def __add__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        self.__support = np.concatenate((self.__support, other.__support))
        self.__values = np.concatenate((self.__values, other.__values))
        self.__aggregate()

    def __aggregate(self):
        """Aggregate values with the same support."""
        self.__support, indices = np.unique(self.__support, return_inverse=True)
        self.__values = np.bincount(indices, weights=self.__values)

    def round_support(self, precision):
        """Round the support of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the support values get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.

        """
        self.__support = np.around(self.__support, precision)
        self.__aggregate()

    def trim(self, cut_off):
        """Trim values below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        self.__support = self.__support[self.__values >= cut_off]
        self.__values = self.__values[self.__values >= cut_off]

    def split_measure(self, cut_off):
        """Split measure into two according to the cut off on values.

        Retain the measure with values greater or equal to the cut off.
        Parameters
        ----------
        cut_off : float
        Returns
        ----------
        other : Measure
            A measure with values strictly below the cut off.

        """
        other = Measure(support=self.__support[self.__values < cut_off],
                        values=self.__values[self.__values < cut_off])
        self.trim(cut_off)
        return other

    def get_cut_off_value(self, P=.99):
        """Get the cut off resulting in optimal P-set.

        Parameters
        ----------
        P : float
            The percentage of the initial total value that the new measure
            will contain. The new measure contains only support with
            heighest values.
        """
        assert 0.0 <= P <= 1.0, "Wrong P for P-optimal set."
        total_value = self.__values.sum()
        i = 0
        S = 0.0
        for intensity in sorted(self.__values):
            S += intensity/total_value
            if S < P:
                break
        return intensity

    def __repr__(self):
        """Represent the measure."""
        pass


class ExperimentalSpectrum(Measure):
    """Store an experimental spectrum."""
    @property
    def mz(self):
        return self.__support

    @property
    def intensity(self):
        return self.__values


class IsotopicDistribution(Measure):
    """Store an isotopic distribution."""
    @property
    def mz(self):
        return self.__support

    @property
    def probability(self):
        return self.__values
