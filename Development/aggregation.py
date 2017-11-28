from math import fsum
import numpy as np
from six.moves import zip


# This is used by the IsotopeCalculator.
# get rid of this when IsoSpec 2.0 is in place
def aggregate_envelopes(masses, probs, digits=2):
    """Aggregate theoretical envelopes.

    Parameters
    ----------
    masses : array
        An array of isotopologues' masses.
    probs : array
        An array of isotopologues' probabilities.
    digits : int
        The number of significant digits used
        while rounding the masses of isotopologues.
    Returns
    ----------
    out : tuple
        A theoretical spectrum of a given resolution.

    """
    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array([k for k in lists])
    newProbs = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs, op_flags=['readwrite']),
                          newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs
