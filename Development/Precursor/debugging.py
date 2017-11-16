%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Data.get_amino_acids import get_amino_acids
from MassTodonPy.MoleculeMaker.Precursor import Precursor

mol = get_dataset('substanceP').precursor
mol.name
mol.fasta
mol.q
mol.modifications

P = Precursor(mol.name,
              mol.fasta,
              mol.q,
              modifications=mol.modifications)

P

get_amino_acids()
P.modifications2
