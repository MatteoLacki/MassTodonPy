%load_ext autoreload
%autoreload
from MassTodon import MassTodon
import numpy as np

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; P = .999

massTodon = MassTodon(fasta, Q)
masses, intensities, noise_masses, noise_intensities = massTodon.randomSpectrum(ionsNo)
mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)
massSpectrum = list(zip(mz,ints))

massTodon.setMassSpectrum(massSpectrum)
