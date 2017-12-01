%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper
from MassTodonPy.MatchMaker.match_cz_ions import match_cz_ions, czMatchMaker, czMatchMakerBasic


# %%time
mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())
D = list(deconvolve(mols, mol.spectrum))

results = [d.report() for d in D]
res = results[0]
M = res['alphas'][0]

matches = czMatchMakerBasic(results, mol.precursor)


options = {'node_color': 'black',
           'node_size': 2,
           'width': 1,
           'with_labels': True,
           'font_size': 10}
nx.draw(matches.graph)
plt.show()
matches.match()
matches.graph.edges



# match_cz_ions(results,
#               mol.precursor.q,
#               mol.precursor.fasta,
#               {},
#               0.0,
#               analyzer='intermediate')
