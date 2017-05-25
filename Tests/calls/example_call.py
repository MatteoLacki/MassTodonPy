from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin

mol = substanceP.copy()
jP  = .999
cutOff = 100
mzPrec = .05
topPercent = .999
precDigits = 2
M_minProb  = .7
max_times_solve = 30
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver  = 'sequential'
method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursorCharge = mol['Q'],
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = mol['modifications']  )
M.readSpectrum( spectrum        = mol['spectrum'],
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent  )

M.prepare_problems(M_minProb)
M.run(  solver  = 'sequential',
        method  = 'MSE',
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
        verbose = verbose )

print M.summarize_results()
print M.gen_ETDetective_inputs()
print M.analyze_reactions('basic')
print M.analyze_reactions('inter')
print M.analyze_reactions('up_inter')
