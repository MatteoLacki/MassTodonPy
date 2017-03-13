%load_ext autoreload
%autoreload
%load_ext line_profiler

from MassTodon import MassTodon
import numpy as np
from pandas import DataFrame

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
fasta='RPKPQQFFGLM'
Q = 3; modifications = {}; ionsNo  = 10000; P = .999

massTodon = MassTodon(fasta, Q, massPrecDigits = 1)
masses, intensities, noise_masses, noise_intensities = massTodon.randomSpectrum(ionsNo)

mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)
massSpectrum = list(zip(mz,ints))

massTodon.peakPicker.setMassSpectrum(massSpectrum)

DF = DataFrame(massTodon.formulator.makeMolecules())
DF.columns = ['type','form','aaNo','q','g']
DF.query("type == 'c71'")

DF.to_csv( path_or_buf="/Users/matteo/Documents/MassTodon/MassTodonPy/Old/compMassTodonR/ubiquitin_formulas_py.csv",
    index = False )
# It seems it's the same.
# DF.query("type == 'c10'")
