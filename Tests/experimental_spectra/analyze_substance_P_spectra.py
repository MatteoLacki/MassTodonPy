import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
import  cPickle as pickle
from    MassTodonPy.MatchMaker import reaction_analist_basic as RA_base, reaction_analist_intermediate as RA_inter, reaction_analist_upper_intermediate as RA_advinter
from    collections import Counter



result_path_specific = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/parsed_sub_P.matteo'
analysis_results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/analyzed_results.json'

with open(result_path_specific,'r') as fp:
    experiments = pickle.load(fp)

def getResults(fasta, Q, WH, WV, L, modifications, spectrum,
               jP=.99, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 100, topPercent = .999, max_times_solve = 10,
               L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False ):
    params = (fasta, Q, WH, WV)

    M = MassTodon(  fasta           = fasta,
                    precursorCharge = Q,
                    precDigits      = precDigits,
                    jointProbability= jP,
                    mzPrec          = mzPrec  )

    M.readSpectrum( spectrum        = spectrum,
                    cutOff          = cutOff,
                    digits          = precDigits,
                    topPercent      = topPercent  )

    M.prepare_problems(M_minProb)

    res = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )
    return params, res

Z = 0.001

%%time
results = [ getResults(*exp,L1_x=Z, L2_x=Z, L1_alpha=Z, L2_alpha=Z) for exp in experiments ]


"_".join(map(str,(fasta, Q, WH, WV)))

analyzers = (RA_base, RA_inter, RA_advinter)
analyzed_results = dict( ("_".join(map(str,(fasta, Q, WH, WV))), map(lambda f:f(res, Q, fasta), analyzers )) for
    (fasta, Q, WH, WV), res in results)

with open(analysis_results_path, 'w') as f:
    json.dump(analyzed_results, f)
