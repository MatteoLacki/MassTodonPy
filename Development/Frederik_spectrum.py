from pandas import DataFrame as df
from MassTodonPy.Data.get_dataset import get_dataset

mol = get_dataset("substanceP")

df({'mz': mol.spectrum.mz, 'intensity':mol.spectrum.intensity}).to_csv(
    path_or_buf='/Users/matteo/Downloads/fwdwsppraca/Frederik_spectrum.csv',
    index=False)
