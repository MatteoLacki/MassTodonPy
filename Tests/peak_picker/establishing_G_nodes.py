from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

mol = substanceP.copy()

cutOff = 100; topPercent = .999; max_times_solve=30
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = True

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

Results = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )

SFG = Results[0]['SFG']
GS = [G for G in SFG.nodes(data=True) if G[1]['type']=='G']

Results[0]
