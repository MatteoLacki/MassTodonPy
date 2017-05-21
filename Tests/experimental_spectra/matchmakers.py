from   MassTodonPy import MassTodon
from   MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from   math import log, lgamma, log10, exp, sqrt
import cPickle as pickle

file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.mzXML'
mass_res_file = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.matteo'

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.99; mzPrec=.05; precDigits=2; M_minProb=.7; cutOff = 100; topPercent = .999
L1_x = L2_x = L1_alpha = L2_alpha = 0.001

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

# modifications = {   'N2':       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }

M.readSpectrum(path=file_path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
M.prepare_problems(M_minProb)

%%time
res = M.run(solver = 'sequential',
            method = 'MSE',
            max_times_solve = 10,
            L1_x   = L1_x, L2_x = L2_x,
            L1_alpha = L1_alpha, L2_alpha = L2_alpha)

with open(mass_res_file, 'w') as f:
    pickle.dump(res, f)

# with open(mass_res_file, 'rb') as f:
#     res = pickle.load(f)

%%time
Prob, Counts = reaction_analist_upper_intermediate(res, Q, fasta, L2=0.0001)
print Prob['ETnoD']
print Prob['ETnoD_prec']
%%time
Prob, Counts = reaction_analist_upper_intermediate(res, Q, fasta, L2=0.0, L1 = 0.1)
print Prob['ETnoD']
print Prob['ETnoD_prec']

%%time
Prob, Counts = reaction_analist_upper_intermediate(res, Q, fasta, L2=0.0, L1 = 0.0)
print Prob['ETnoD']
print Prob['ETnoD_prec']


%%time
Prob, Counts = reaction_analist_intermediate(res, Q, fasta)
print Prob['ETnoD']
print Prob['ETnoD_prec']


%%time
EasyRes = reaction_analist_basic(res, Q, fasta)
EasyRes['PTR']
EasyRes['ETnoD']
