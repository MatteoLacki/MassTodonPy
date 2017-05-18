from    MassTodonPy  import MassTodon
from    collections  import Counter
from    time         import time

jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
cutOff = 100.; cutOff2=0.0; topPercent = .999; max_times_solve=30
L1_x=0.001; L2_x=0.001; L1_alpha=0.001; L2_alpha=0.001; verbose=False

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8
file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.mzXML'

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec  )

M.readSpectrum( path        = file_path,
                cutOff      = cutOff,
                digits      = precDigits,
                topPercent  = topPercent    )

M.prepare_problems(M_minProb)

T0_deconv = time()
Results = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )
T1_deconv = time()
T_deconv  = T1_deconv - T0_deconv


M.analyze_reactions(analyzer='basic')
M.analyze_reactions(analyzer='inter')
M.analyze_reactions(analyzer='up_inter')
