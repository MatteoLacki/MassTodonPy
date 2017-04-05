%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import  MassTodon
import  cPickle         as      pickle
# from    Formulator      import  makeFormulas
# import  networkx        as      nx

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# modifications = {   ('N',2) :       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }
# Forms = makeFormulas(fasta=fasta, Q=Q, fragType='cz')
# M = MassTodon(  fasta           = fasta,
#                 precursorCharge = Q,
#                 precDigits      = precDigits,
#                 jointProbability= jP,
#                 mzPrec          = mzPrec )
# path  = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
# path += 'Ubiquitin_ETD_10 ms_1071.mzXML'
# cutOff = 100; topPercent = .999
# M.readSpectrum(path=path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
# M.prepare_problems(M_minProb)
# mu=1e-5; lam=0.0; nu=0.001
# %%time
# res = M.run(solver='sequential', method='MSE', mu=mu, lam=lam, nu=0.001)
# with open(file_path, 'w') as f:
    # pickle.dump(res, f)
with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

from MatchMaker import reaction_analist_basic

%%time
results = reaction_analist_basic(MassTodonResults, fasta, Q) # works perfectly!!!
results
