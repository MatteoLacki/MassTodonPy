%load_ext autoreload
%autoreload 2

import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch

mol = get_dataset('substanceP')
mols = list(mol.precursor.molecules())
str(mols[0].formula)

D = list(deconvolve(mols, mol.spectrum))
results = [d.report() for d in D]
res = results[0]
M = res['alphas'][0]

# matches = SimpleCzMatch(results, mol.precursor)
# matches.get_probabilities()
# matches.get_intensities()



nx.max_flow_min_cost
