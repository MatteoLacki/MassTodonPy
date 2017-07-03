from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts import substanceP, ubiquitin
from time  import time

mol = substanceP.copy()
# mol = ubiquitin.copy()

jP      = .999
mzPrec  = .05
opt_P   = .99
max_times_solve = 10
multiprocesses_No = None
verbose = True
solver  = 'multiprocessing'
# solver  = 'sequential'

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = mzPrec,
                    joint_probability_of_envelope= jP,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = opt_P,
                    solver          = solver,
                    multiprocesses_No = multiprocesses_No,
                    max_times_solve = max_times_solve,
                    verbose         = verbose )


print res
