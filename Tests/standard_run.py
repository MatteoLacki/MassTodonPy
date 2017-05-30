from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts import substanceP, ubiquitin
import pandas as pd


mol = substanceP.copy()
Results = MassTodonize(
    fasta           = mol['fasta'],
    precursor_charge= mol['Q'],
    joint_probability_of_envelope = .990,
    mz_prec         = .05,
    modifications   = mol['modifications'],
    spectrum        = mol['spectrum'],
    opt_P           = 1.00,
    min_prob_of_envelope_in_picking = .7,
    solver          = 'sequential',
    method          = 'MSE',
    max_times_solve = 10,
    L1_x            = .001,
    L2_x            = .001,
    L1_alpha        = .001,
    L2_alpha        = .001,
    verbose         = False )
