from MassTodonPy import MassTodon, MassTodonize
from MassTodonPy.TestScripts import substanceP, ubiquitin
from time  import time

# mol = substanceP.copy()
mol = ubiquitin.copy()

jP = .999

cutOff = 100
mzPrec = .05
opt_P  = .99
M_minProb  = .7
max_times_solve = 10

# solver  = 'sequential'
solver  = 'multiprocessing'
multiprocesses_No = 4

method  = 'MSE'
verbose = True


T0 = time()

#
# M = MassTodon(  fasta           = mol['fasta'],
#                 precursor_charge= mol['Q'],
#                 joint_probability_of_envelope= jP,
#                 mz_prec         = mzPrec,
#                 modifications   = mol['modifications'],
#                 verbose         = verbose   )
# #
# M.read_n_preprocess_spectrum(
#     spectrum = mol['spectrum'],
#     cut_off  = cutOff  )
#
# M.run(  solver  = solver,
#         multiprocesses_No = multiprocesses_No,
#         min_prob_per_molecule = M_minProb,
#         method  = method,
#         max_times_solve = max_times_solve )

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = mzPrec,
                    joint_probability_of_envelope= jP,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = opt_P,
                    solver          = solver,
                    multiprocesses_No = 4,
                    max_times_solve = max_times_solve,
                    verbose         = verbose )


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
