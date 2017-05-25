import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon

# fasta, Q, WH, WV, L, modifications, spectrum = exp
# cutOff = 100; topPercent = .999; max_times_solve=30
# jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# L1_x = L2_x = L1_alpha = L2_alpha = .001
# verbose = True
def getResults(fasta, Q, WH, WV, L, modifications, spectrum, jP=.999, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 100., cutOff2=0.0, topPercent = .999, max_times_solve=30, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False):
    params = (fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha)
    try:
        M = MassTodon(  fasta           = fasta,
                        precursorCharge = Q,
                        precDigits      = precDigits,
                        jointProbability= jP,
                        mzPrec          = mzPrec,
                        modifications   = modifications  )

        M.readSpectrum( spectrum        = spectrum,
                        cutOff          = cutOff,
                        digits          = precDigits,
                        topPercent      = topPercent  )

        M.prepare_problems(M_minProb)
        T0_deconv = time()
        Results   = M.run(solver  = 'sequential',
                        method  = 'MSE',
                        max_times_solve = max_times_solve,
                        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                        verbose = verbose )
        T1_deconv = time()
        T_deconv  = T1_deconv - T0_deconv
        RA = {}

        try:
            RA['base'] = M.analyze_reactions(analyzer='basic')
        except:
            print 'Base missing'
        try:
            RA['inter'] = M.analyze_reactions(analyzer='inter')
        except:
            print 'Intermediate missing'
        try:
            RA['up_inter'] = M.analyze_reactions(analyzer='up_inter')
        except:
            print 'Upper Intermediate missing'

        res = (params, Results, WH, WV, RA)
    except Exception as e:
        res = (params, e)
    return res
