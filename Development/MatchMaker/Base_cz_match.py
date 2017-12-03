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
matches.lavish
matches.ETD
matches.broken_bond
matches.ETnoD
matches.PTR




def __get_probs(counts, Probs, tag1, tag2, name1=None, name2=None):
    """Make two probabilities out of counts and name them properly."""
    if not name1:
        name1 = tag1
    if not name2:
        name2 = tag2
    total = counts[tag1]+counts[tag2]
    if total > 0.0:
        Probs[name1] = float(counts[tag1])/total
        Probs[name2] = 1.0 - Probs[name1]
    return Probs
