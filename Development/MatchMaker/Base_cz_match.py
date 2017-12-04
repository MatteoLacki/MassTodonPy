%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper
from MassTodonPy.MatchMaker.Base_cz_match import Base_cz_match


# %%time
mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())
D = list(deconvolve(mols, mol.spectrum))

results = [d.report() for d in D]
res = results[0]
M = res['alphas'][0]

matches = Base_cz_match(results, mol.precursor)
matches.get_probabilities()
x = matches.get_intensities()
