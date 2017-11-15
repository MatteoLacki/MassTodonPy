import pkg_resources
import json
from collections import defaultdict


def get_elements():
    return {
        "O", "Xe", "Cs", "Hg", "S", "Ru", "H", "Zn", "Sr", "Al", "Sm",
        "Zr", "Ho", "Ta", "Pb", "Te", "He", "Ti", "As", "Ge", "Pr", "U",
        "Tl", "Ir", "Tm", "Fe", "Si", "Cl", "Eu", "Tb", "W", "Er", "P",
        "Os", "K", "Dy", "Lu", "Bi", "Ga", "Pt", "La", "Be", "F", "Yb",
        "Kr", "Cd", "Mn", "Ar", "Cr", "Se", "Sb", "Hf", "Sc", "Ca", "Ba",
        "Rb", "Sn", "Co", "Cu", "Ne", "Pd", "In", "N", "Au", "Y", "Ni",
        "Rh", "C", "Li", "Th", "B", "Mg", "Na", "Pa", "V", "Re", "Nd",
        "Br", "Ce", "I", "Ag", "Gd", "Nb", "Mo"}


def get_isotopic_masses_and_probabilities_raw():
    """
    Retrieve the information on masses and frequencies of isotopes.

    Returns
    =======
    res : list
    """
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path + "isotopes.json", "rb") as f:
        isotopes = json.load(f)
    return isotopes


def get_isotopic_masses_and_probabilities():
    """
    Retrieve the information on masses and frequencies of isotopes.

    Returns
    =======
    Two dictionaries, both with keys set to element encodings,
    such as Ag, Al, Ar, ...
    The values of the first dictionary are lists of masses of elements.
    The values of the first dictionary are lists of probabilities of elements.
    """
    isotopes_raw = get_isotopic_masses_and_probabilities_raw()
    iso_masses = defaultdict(list)
    iso_probs = defaultdict(list)
    for element, isos in isotopes_raw:
        for mass, prob in isos:
            iso_masses[element].append(mass)
            iso_probs[element].append(prob)
    return iso_masses, iso_probs
