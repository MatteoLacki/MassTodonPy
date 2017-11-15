%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_data import get_dataset, get_amino_acids
from MassTodonPy.MoleculeMaker.Precursor import Precursor

mol = get_dataset('substanceP').precursor

P = Precursor(mol.name,
              mol.fasta,
              mol.q,
              modifications=mol.modifications)

get_amino_acids()
P.modifications2
