%load_ext autoreload
%autoreload 2

from math import ceil, log10
import networkx as nx
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.PeakPicker.PeakPicker2 import get_deconvolution_problems

mol = get_dataset('substanceP')

%%timeit
molecules = mol.precursor.molecules()
problems = list(get_deconvolution_problems(molecules, mol.spectrum))
