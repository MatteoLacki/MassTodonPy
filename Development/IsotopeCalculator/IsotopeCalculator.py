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
z11.__dict__

str(z11.formula)


subP.spectrum.round_masses(1)

list(subP.spectrum[300,1000])

.reset_parser('([a-z]?)([0-9])')

str(z11.formula)
str(z11.formula)


%%timeit
iso_calc = IsotopeCalculator(mz_precision=2)
spectrum = iso_calc.get_envelope(z11.formula, .99, z11.q, z11.g)

iso_calc.get_mean_mz(z11.formula, z11.q, z11.g)
iso_calc.get_monoisotopic_mz(z11.formula, z11.q, z11.g)
iso_calc.get_mz_sd(z11.formula, z11.q, z11.g)
