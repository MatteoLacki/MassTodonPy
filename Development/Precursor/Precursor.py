%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Precursor.Precursor import Precursor

subP = get_dataset('substanceP')

precursor = subP.precursor
MM = precursor.molecules()
mols = list(MM)
