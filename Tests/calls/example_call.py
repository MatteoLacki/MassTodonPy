# %load_ext autoreload
# %autoreload
import os

# Setting up locally the number of threads used by
# BLAS to 1.
old = os.environ.get('OMP_NUM_THREADS', None)
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
from    time  import time

# mol = substanceP.copy()
mol = ubiquitin.copy()

jP  = .999
cutOff = 100
mzPrec = .05
M_minProb  = .7
max_times_solve = 10
L1_x = L2_x = L1_alpha = L2_alpha = .001
# solver  = 'sequential'
solver  = 'multiprocessing'
method  = 'MSE'
verbose = True


T0 = time()
M = MassTodon(  fasta           = mol['fasta'],
                precursor_charge= mol['Q'],
                joint_probability_of_envelope= jP,
                mz_prec         = mzPrec,
                modifications   = mol['modifications'],
                verbose         = verbose   )

M.read_n_preprocess_spectrum(
    spectrum = mol['spectrum'],
    cut_off  = cutOff  )

M.prepare_problems(M_minProb)

M.run(  solver  = solver,
        method  = method,
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha)
T1 = time()

print 'Envelopes Generation', M.IsoCalc.stats['Envelopes Generation Total T']
print
print 'Total Time', T1-T0
print

# print M.summarize_results()
# print M.gen_ETDetective_inputs()
# print M.analyze_reactions('basic')
# print M.analyze_reactions('intermediate')
# print M.analyze_reactions('advanced')

if old:
    os.environ['OMP_NUM_THREADS'] = old
else:
    del os.environ['OMP_NUM_THREADS']
