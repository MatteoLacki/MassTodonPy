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

def getResults(fasta, Q, WH, WV, L, modifications, spectrum,
               jP=.99, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 0, topPercent = .999,
               L1_x=0.0, L2_x=0.0, L1_alpha=0.0, L2_alpha=0.0 ):
    params = (fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, L1_x, L2_x, L1_alpha, L2_alpha)
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
                    L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                    verbose = True )
    return res


results = [ getResults(*exp) for exp in experiments ]
res2    = [ (alphas, error, sol, params, SFG) for r in results for alphas, error, sol, params, SFG in r if sol['status'] != 'optimal']

# which spectrum has the most?
all_sols=[]
for r in results:
    sols = Counter()
    for alphas, error, sol, params, SFG in r:
        sols[sol['status']]+=1
    all_sols.append(sols)

[ n for n, sol in enumerate(all_sols) if sol['unknown']==3 ]
# spectra 5, 25, 28



fasta, Q, WH, WV, L, modifications, spectrum = experiments[25]
L1_x=0.0; L2_x=0.0; L1_alpha=0.0; L2_alpha=0.0

res = getResults(   fasta, Q, WH, WV, L, modifications, spectrum,
    L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha )

alphas, error, sol, params, SFG = res[0]
P, q, G, h, A, b, x0 = params
from cvxopt import matrix, spmatrix, sparse, spdiag, solvers

def kkt(W):
    print W

solvers.qp(P, q, G, h, A, b, initvals=x0)

solvers.qp(P, q, G, h, A, b, initvals=x0, kktsolver = kkt)

%%time
solvers.qp(P, q, G, h, A, b, initvals=x0, kktsolver = 'ldl2')

%%time
solvers.qp(P, q, G, h, A, b, initvals=x0, kktsolver = 'chol')


print A
print A[4,1]




SFG.nodes(data=True)


Counter(sol['status'] for r in res)

alphas, error, sol, params, SFG =  res[1][0]


sum([len(r[1]) for r in results])
len(res2)

res2[0]

# res2[0]
# alphas, error, sol, params, SFG = res2[0]
#
# P, q, G, h, A, b = params
# print P
# print q
# print q
