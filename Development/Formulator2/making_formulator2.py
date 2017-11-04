%load_ext autoreload
%autoreload 2

from MassTodonPy import get_dataset
from MassTodonPy.Formulator.formulator import  make_formulas
mol = get_dataset('substanceP')

Forms = make_formulas(
    fasta=mol["fasta"],
    Q=mol["Q"],
    frag_type="cz",
    distance_charges=5,
    modifications=mol["modifications"])

list(Forms.makeMolecules())
