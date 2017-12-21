%load_ext autoreload
%autoreload 2

from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolution.Deconvolve import deconvolve

# %time
mol = get_dataset('substanceP') # adjust the spectrum
molecules = list(mol.precursor.molecules())
sigma = 0.01949749
support_length = 0.1

deconvolution_graph = deconvolve(molecules,
                                 mol.spectrum,
                                 'Ciacho_Wanda',
                                 mz_tol=.05,
                                 min_prob_per_molecule=.7)

problems = list(deconvolution_graph)
problem = problems[0]
problem.node['M0']
