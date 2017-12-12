%load_ext autoreload
%autoreload 2

from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.PeakPicker import get_deconvolution_problems
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve

# %%time
mol = get_dataset('substanceP') # adjust the spectrum
molecules = list(mol.precursor.molecules())
sigma = 0.01949749
support_length = 0.1

deconvolution_graph = deconvolve(molecules,
                                 mol.spectrum,
                                 mz_tol=support_length/2,
                                 method='Ciacho_Wanda')

problems = list(deconvolution_graph)

problems[0].nodes


def fun(arg=0, **args):
    print(arg, args)

fun(**{'a': 10, 'b': 20})

deconvolve()

problems[0].nodes()
nx.draw(problems[0], node_labels=True)
plt.show()

Counter(n[0] for n in problem[0].nodes())
