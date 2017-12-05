%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch
from MassTodonPy.MatchMaker.CzMatch import CzMatch
from MassTodonPy.MatchMaker.RegularizedCzMatch import RegularizedCzMatch

# %%time
mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())

strings = []
for m in mols:
    if m.name is not 'precursor':
        strings.append("({}_{{{}}}, {}), ".format(m.name[0], m.name[1:], m.q))

"".join(set(strings))

strings_full = []
for m in mols:
    if m.name is not 'precursor':
        strings_full.append("({}_{{{}}}, {}, {}), ".format(m.name[0],
                                                           m.name[1:],
                                                           m.q, m.g))

"".join(set(strings_full))
