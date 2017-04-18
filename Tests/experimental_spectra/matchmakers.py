%load_ext autoreload
%load_ext line_profiler
%autoreload


from   MassTodonPy import MassTodon
from   MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from   math import log, lgamma, log10, exp, sqrt
import cPickle as pickle

file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.mzXML'

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7; cutOff = 100; topPercent = .999
mu=1e-5; lam=0.0; nu=0.001

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

# modifications = {   ('N',2) :       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }

M.readSpectrum(path=file_path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
M.prepare_problems(M_minProb)


# %%time
res = M.run(solver='sequential', method='MSE', mu=mu, lam=lam, nu=0.001)
MassTodonResults = res
# with open(file_path, 'w') as f:
    # pickle.dump(res, f)
# with open(file_path, 'rb') as f:
#     MassTodonResults = pickle.load(f)

%%time
LogProb = reaction_analist_advanced(MassTodonResults, Q, fasta, maxIter=100, const=1000, eps = 0.0, verbose=True)


%%time
EasyRes = reaction_analist_basic(MassTodonResults, Q, fasta)

EasyRes['PTR']
EasyRes['ETnoD']
