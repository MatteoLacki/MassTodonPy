from    MassTodonPy  import MassTodon
from    collections  import Counter, defaultdict
from    read_experiments import experiments as substancesP

jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
cutOff = 100.; cutOff2=0.0; topPercent = .999; max_times_solve=30
L1_x=0.001; L2_x=0.001; L1_alpha=0.001; L2_alpha=0.001; verbose=False


X = defaultdict(list)
for i, (_,_,WH,WV,_,modifications,_) in enumerate(substancesP):
    X[(WH,WV)].append((i,modifications))


specNo = X[(150,700)][0][0]
specNo = X[(125,300)][0][0]
fasta, Q, WH, WV, L, modifications, spectrum = substancesP[specNo]

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = modifications )

M.readSpectrum( spectrum        = spectrum,
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent    )

M.prepare_problems(M_minProb)

Results = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )


M.analyze_reactions(analyzer='basic')
M.analyze_reactions(analyzer='inter')
M.analyze_reactions(analyzer='up_inter')
