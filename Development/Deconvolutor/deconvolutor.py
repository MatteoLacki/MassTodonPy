%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper

%%time
mol = get_dataset('substanceP') # adjust the spectrum
mols = list(mol.precursor.molecules())
D = list(deconvolve(mols, mol.spectrum))


# TODO play with multiprocessing later on.
# %%time
# mol = get_dataset('substanceP') # adjust the spectrum
# mols = list(mol.precursor.molecules())
# D, runtime = deconvolve(mols, mol.spectrum, processes=3)
