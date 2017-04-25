import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
import  cPickle as pickle
from    MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from    collections import Counter

result_path_specific = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/parsed_sub_P.matteo'

with open(result_path_specific,'r') as fp:
    experiments = pickle.load(fp)

def getResults(fasta, Q, WH, WV, L, modifications, spectrum, jP=.999, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 100, topPercent = .999, mu=1e-5, lam=0.0, nu=0.001):
    params = (fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, mu, lam, nu)
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
    Results = M.run(solver  = 'sequential',
                    method  = 'MSE',
                    mu=mu, lam=lam, nu=0.001,
                    verbose = True )
    res = (True, Results)
    return res

%time
results = [ getResults(*exp) for exp in experiments ]

Counter(r[0] for r in results)
res2 = [ (alphas, error, sol, params, SFG) for r in results if r[0] for alphas, error, sol, params, SFG in r[1] if sol['status'] != 'optimal']

sum([len(r[1]) for r in results])
len(res2)

res2
Stati = Counter( sol['status'] for alphas, error, sol, params, SFG in res2)
Stati
