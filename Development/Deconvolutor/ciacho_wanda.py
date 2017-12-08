%load_ext autoreload
%autoreload 2

from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.PeakPicker import get_deconvolution_graphs


# %%time
mol = get_dataset('substanceP') # adjust the spectrum


molecules = list(mol.precursor.molecules())
support_length = .2

deconvolution_graph = get_deconvolution_graphs(molecules,
                                               mol.spectrum,
                                               mz_tol=support_length/2)

problem = list(deconvolution_graph)



problem[0].nodes()
nx.draw(problem[0], node_labels=True)
plt.show()

Counter(n[0] for n in problem[0].nodes())
