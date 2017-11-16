%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_amino_acids import get_amino_acids
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.MoleculeMaker import MoleculeMaker

subP = get_dataset('substanceP')
precursor = subP.precursor

MM = MoleculeMaker(precursor)

list(MM.uncharged_molecules())
list(MM.charged_molecules())
