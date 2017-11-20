%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.MoleculeMaker import MoleculeMaker, Molecule

subP = get_dataset('substanceP')
precursor = subP.precursor

MM = MoleculeMaker(precursor)

list(MM.charged_molecules())
