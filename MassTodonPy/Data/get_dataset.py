import pkg_resources
import json
import numpy as np
from collections import namedtuple

from MassTodonPy.MoleculeMaker.Precursor import Precursor
from MassTodonPy.Spectra.Spectrum import Spectrum


Dataset = namedtuple('Dataset', 'precursor spectrum instrument')


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

    spectrum = Spectrum(*tuple(np.array(d) for d in mol['spectrum']))

    modifications = {int(k): v for k, v in
                     mol['modifications'].items()}

    precursor = Precursor(name=mol['name'],
                          fasta=mol['fasta'],
                          q=mol['Q'],
                          modifications=modifications)

    instrument = {}
    if dataset_name is 'substanceP':
        instrument['name'] = 'synapt'
        instrument['wave height'] = 0
        instrument['wave velocity'] = 300

    elif dataset_name is 'ubiquitin':
        instrument['name'] = 'orbitrap'
        instrument['acquisition time'] = '10 ms'

    return Dataset(precursor=precursor,
                   spectrum=spectrum,
                   instrument=instrument)
