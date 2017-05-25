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

fasta,Q,WH,WV,L,modifications,spectrum = experiments[25]


jP=.99; mzPrec=.05; precDigits=2; M_minProb=.7; cutOff = 0
topPercent = .999;
L1_x=0.0; L2_x=0.0; L1_alpha=0.0; L2_alpha=0.0

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
res = M.run(    solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = 5,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = False )

res

reaction_analist_upper_intermediate(res, Q, fasta, verbose=False)
reaction_analist_intermediate(res, Q, fasta, verbose=False)
reaction_analist_basic(res, Q, fasta)

res[0]['alphas']

Counter()
res

reaction_analist_upper_intermediate
