%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import Deconvolutor

mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())
D = Deconvolutor(mols, mol.spectrum)

for d in D.problems:
    print(d.report()['status'])
