%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_data import get_isotopic_masses_and_probabilities
import json
from collections import namedtuple, defaultdict

data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"

masses, probs = get_isotopic_masses_and_probabilities()

isotopes = [(aa, list(zip(masses[aa], probs[aa]))) for aa in masses]

with open(data_path + "isotopes.json", "w") as f:
    json.dump(isotopes, f)
