%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.MoleculeMaker import get_molecules
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator

subP = get_dataset('substanceP')
precursor = subP.precursor
mols = get_molecules(precursors=[precursor])

for mol in mols:
    print(mol)

z11 = mol

iso_calc = IsotopeCalculator()
spectrum = iso_calc.get_envelope(.99, z11.formula, z11.q, z11.g)


spectrum.mz
spectrum.probability
