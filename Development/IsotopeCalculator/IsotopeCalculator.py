%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.MoleculeMaker import molecules
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator

subP = get_dataset('substanceP')
precursor = subP.precursor
mols = molecules(precursors=[precursor])

subP.spectrum.mass
subP.spectrum.intensity

for mol in mols:
    print(mol)

prec = mol

iso_calc = IsotopeCalculator()
spectrum = iso_calc.get_envelope(.99, prec.formula, prec.q, prec.g)


spectrum.mass
spectrum.probability

x = ""
assert x is not ""
