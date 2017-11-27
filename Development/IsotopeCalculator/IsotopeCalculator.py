%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator

subP = get_dataset('substanceP')
precursor = subP.precursor
mols = precursor.molecules()
for mol in mols:
    mol
z11 = mol

str(z11.formula)


%%timeit
iso_calc = IsotopeCalculator(mz_precision=3)
spectrum = iso_calc.get_envelope(z11.formula, .99, z11.q, z11.g)
