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

mol.precursor.formula

mols = list(mol.precursor.molecules())
str(mols[0].formula)

D = list(deconvolve(mols, mol.spectrum))
results = [d.report() for d in D]
res = results[0]
M = res['alphas'][0]

# matches = SimpleCzMatch(results, mol.precursor)
# matches.get_probabilities()
# matches.get_intensities()

matches = CzMatch(results, mol.precursor)
matches.get_probabilities()
matches.get_intensities()

matches.graph.edges.data('flow')
