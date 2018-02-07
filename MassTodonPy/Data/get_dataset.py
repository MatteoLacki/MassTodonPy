import pkg_resources
import json
import numpy as np

from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.Spectra.Spectrum import Spectrum


class Dataset(object):
    def __init__(self, precursor, spectrum, instrument):
        self.precursor = precursor
        self.spectrum = spectrum
        self.instrument = instrument

    def __repr__(self):
        out = "---- Dataset ----\n{}\n".format(self.precursor.__repr__())
        out += "{}Instrument{}\n".format(self.spectrum.__repr__(),
                                         self.instrument.__repr__())
        out += "----------------"
        return out

def get_dataset(dataset_name, json=False):
    """
    Retrieve examplary spectra informations.

    Parameters
    ==========
    dataset_name: string
        Warning
        =======
        Can take values: substanceP or ubiquitin.
    json : boolean
        Should we read in json-saved file instead of python dictionary?

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
    assert dataset_name in ['substanceP',]
    # ['substanceP', 'ubiquitin'] # TODO fix ubiquitin example

    if json:
        path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
        with open(path + dataset_name + '.json', 'rb') as f:
            mol = json.load(f)
    else:
        if dataset_name is 'substanceP':
            from MassTodonPy.Data.substanceP import substanceP as mol
        else:
            raise AttributeError("No data set of that name.")
    spectrum = Spectrum(spectrum=mol['spectrum'])

    modifications = {int(k): v for k, v in
                     mol['modifications'].items()}

    precursor = Precursor(name=mol['name'],
                          fasta=mol['fasta'],
                          charge=mol['Q'],
                          modifications=modifications,
                          fragmentation_type="cz")

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
