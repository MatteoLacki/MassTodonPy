%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset

substanceP = get_dataset('substanceP')
substanceP.spectrum.plot()
