%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.MoleculeMaker import MoleculeMaker, Molecule
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator

subP = get_dataset('substanceP')
precursor = subP.precursor

MM = MoleculeMaker(precursor)
molecules = list(MM.charged_molecules())
prec = molecules[5]

iso_calc = IsotopeCalculator()
spectrum = iso_calc.get_envelope(.99, prec.formula, prec.q, prec.g)


spectrum.mass
spectrum.probability
