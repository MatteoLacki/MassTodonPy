import pkg_resources
import json
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle


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
    with open(path + dataset_name + ".json", "rb") as f:
        data = json.load(f)
    data["spectrum"] = tuple(np.array(d) for d in data["spectrum"])
    return data


def get_bricks():
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
        bricks = pickle.load(f)
    return bricks


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
