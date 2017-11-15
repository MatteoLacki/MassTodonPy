%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset

subP = get_dataset('substanceP')

subP.precursor
