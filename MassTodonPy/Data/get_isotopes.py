import pkg_resources
import json
from collections import defaultdict


def get_isotopic_masses_and_probabilities():
    """Retrieve the information on masses and frequencies of isotopes.

    Returns
    =======
    out : tuple
        Two dictionaries, both with keys set to element encodings,
        such as Ag, Al, Ar, ...
        The values of the first dictionary are lists of masses of elements.
        The values of the first dictionary are lists of probabilities of elements.
    """
    from MassTodonPy.Data.isotopes import isotopes as isotopes_raw
    iso_masses = defaultdict(list)
    iso_probs = defaultdict(list)
    for element, isos in isotopes_raw:
        for mass, prob in isos:
            iso_masses[element].append(mass)
            iso_probs[element].append(prob)
    return iso_masses, iso_probs
