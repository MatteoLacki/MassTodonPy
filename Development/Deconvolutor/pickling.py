%load_ext autoreload
%autoreload 2

import pickle
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.PeakPicker import get_deconvolution_graphs
from MassTodonPy.Deconvolutor.DeconvolutionProblem import DeconvolutionProblem

mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())

graphs = list(get_deconvolution_graphs(mols, mol.spectrum))
graph = graphs[0]



DP = DeconvolutionProblem(graph)

DP.solve()
DP.report()


with open('nozyce_dupa_chuj.pickle', 'wb') as f:
    pickle.dump(graph, f)


with open('nozyce_dupa_chuj.pickle', 'wb') as f:
    pickle.dump(DP, f)


list(graph.nodes.items())

for a, b in graph.nodes.data(('mz','intensity')):
    print(a, b)


a, b, c = 1, *[2,3]
