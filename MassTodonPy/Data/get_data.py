import pkg_resources
import json
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle


element_tags = {
    "O", "Xe", "Cs", "Hg", "S", "Ru", "H", "Zn", "Sr", "Al", "Sm",
    "Zr", "Ho", "Ta", "Pb", "Te", "He", "Ti", "As", "Ge", "Pr", "U",
    "Tl", "Ir", "Tm", "Fe", "Si", "Cl", "Eu", "Tb", "W", "Er", "P",
    "Os", "K", "Dy", "Lu", "Bi", "Ga", "Pt", "La", "Be", "F", "Yb",
    "Kr", "Cd", "Mn", "Ar", "Cr", "Se", "Sb", "Hf", "Sc", "Ca", "Ba",
    "Rb", "Sn", "Co", "Cu", "Ne", "Pd", "In", "N", "Au", "Y", "Ni",
    "Rh", "C", "Li", "Th", "B", "Mg", "Na", "Pa", "V", "Re", "Nd",
    "Br", "Ce", "I", "Ag", "Gd", "Nb", "Mo"}


def get_dataset(dataset_name):
    """
    Retrieve examplary spectra informations.

    Parameters
    ==========
    dataset_name: string
        Warning
        =======
        Can take values: substanceP or ubiquitin.

    Returns
    =======
    A dictionary with fields:
    Q: int
        The charge of the precursor.
    WH: int
        Wave height.
    WV: int
        Wave velocity.
    fasta: string
        The chain of amino acids that make up the precursor.
    modifications: dictionary of dictionaries
        Atom diffs of modifications of the atomic content.
    name: string
        The name of the dataset.
    spectrum:
        A list of lists of floats corresponding to masses and intensities.
    """
    assert dataset_name in ['substanceP', 'ubiquitin']
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path + dataset_name + '.json', 'rb') as f:
        mol = json.load(f)
    mol['spectrum'] = tuple(np.array(d) for d in mol['spectrum'])
    mol['modifications'] = {int(k): v for k, v in mol['modifications'].items()}
    return mol


def get_amino_acids():
    """
    Retrieve the information on amino acidic bricks.

    Returns
    =======
    A dictionary with keys set to element encodings,
    such as Ag, Al, Ar, ...
    The values are yet again dictionaries, with keys L, C, R,
    that stand for left, center and right.
    L corresponds to the atoms next to nitrogen,
    C to atoms next to C_alpha,
    and R to atoms next to the carbon from the carboxyl group.
    Values correspond to atom counts gathered in a linear counter.
    """
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path+"amino_acids.pickle", "rb") as f:
        amino_acids = pickle.load(f)
    return amino_acids


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
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path + "isotopes.pickle", "rb") as f:
        iso_masses, iso_probs = pickle.load(f)
    return iso_masses, iso_probs
