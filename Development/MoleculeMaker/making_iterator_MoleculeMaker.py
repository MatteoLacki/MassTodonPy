%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_data import get_dataset, get_amino_acids
from MassTodonPy.MoleculeMaker.MoleculeMaker import MoleculeMaker, make_molecules

mol = get_dataset('substanceP')
Molecules = MoleculeMaker(mol.precursor)

# Something is really wrong with the bloody modification parser.
# Keeps on adding new empty structures


mol.precursor.modifications
list(Molecules)


amino_acids = get_amino_acids()
Molecules.make_superatoms()
L = Molecules.superAtoms.copy()
L.reverse()

amino_acids[mol.precursor.fasta[0]]


NC = Molecules.get_amino_acids(direction="N -> C")
NC_L = list(NC)
Molecules.superAtoms
Molecules.superAtoms == NC_L

CN = Molecules.get_amino_acids(direction="N <- C")
L
CN_L = list(CN)
CN_L

NC_L.reverse()
NC_L == CN_L

mol.precursor.modifications
