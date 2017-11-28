%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator

subP = get_dataset('substanceP')

precursor = subP.precursor
mols = precursor.molecules()
mols = list(mols)
z11 = mols[-1]

z11.monoisotopic_mz
z11.mean_mz
z11.isotopologues
