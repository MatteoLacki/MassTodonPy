from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
from MassTodonPy.Formulator.bricks import makeBricks


mol = substanceP.copy()
jP  = .999
cutOff = 100
mzPrec = .05
topPercent = .999
M_minProb  = .7
max_times_solve = 30
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver  = 'sequential'
method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursor_charge= mol['Q'],
                joint_probability_of_envelope= jP,
                mz_prec         = mzPrec,
                modifications   = mol['modifications']  )

M.read_n_preprocess_spectrum(
    spectrum = mol['spectrum'],
    cut_off  = cutOff  )

import cPickle as pickle

bricks = makeBricks()
bricks_file = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/amino_acids.txt'

with open(bricks_file, 'w') as f:
    pickle.dump(bricks, f)
