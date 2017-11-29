%load_ext autoreload
%autoreload 2

from inspect import getsourcelines
from collections import Counter
from math import ceil, log10
import networkx as nx
from networkx import connected_component_subgraphs, connected_components

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.PeakPicker.PeakPicker2 import get_deconvolution_problems

mol = get_dataset('substanceP')

molecules = mol.precursor.molecules()
DG = get_deconvolution_problems(molecules, mol.spectrum)
Counter(n[0] for n in DG.nodes)

CC = list(connected_component_subgraphs(DG))
len(CC)



%%timeit
molecules = mol.precursor.molecules()
# problems = list(get_deconvolution_problems(molecules, mol.spectrum))






help(connected_component_subgraphs)
help(connected_components)
getsourcelines(connected_component_subgraphs)
