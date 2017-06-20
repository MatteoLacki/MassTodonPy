from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin

# mol = substanceP.copy()
mol = ubiquitin.copy()

jP  = .999
cutOff = 100
mzPrec = .05
topPercent = .999
M_minProb  = .7
max_times_solve = 10
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver  = 'sequential'
method  = 'MSE'
verbose = True

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

M.run(  solver  = 'sequential',
        method  = 'MSE',
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha)

print M.summarize_results()

print M.gen_ETDetective_inputs()
print M.analyze_reactions('basic')
print M.analyze_reactions('inter')
print M.analyze_reactions('up_inter')
