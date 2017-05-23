from run_massTodon_on_WH_WV import getResults
from collections import Counter, defaultdict
import numpy as np
import numpy.random as npr
import cPickle as pickle
from  read_experiments import experiments as substancesP

len(substancesP)
len(substancesP[0])

fasta, Q, WH, WV, L, modifications, spectrum = substancesP[0]
verbose = True

params, Results, WH, WV, RA = getResults(fasta, Q, WH, WV, L, modifications, spectrum, verbose=verbose)
fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha = params


SFG = Results[0]['SFG']

from interval import interval as P, inf

for G in SFG:
    M_G = Counter()
    mz_interval = P[0.0, inf]
    if SFG.node[G]['type']=='G':
        for I in SFG[G]:
            I_mz = SFG.node[I]['mz']
            mz_interval = mz_interval & P[I_mz-mzPrec, I_mz+mzPrec]
            for M in SFG[I]:
                if SFG.node[M]['type']=='M':
                    M_G[ tuple(SFG.node[M][k] for k in ('molType','formula','q','g')) ] += SFG.edge[G][I]['estimate']
        print G, mz_interval





path_iter = GIM_paths(SFG)
G, I, M = next(path_iter)

SFG.node[G]
SFG.edge[G][I]['estimate']
