%load_ext autoreload
%load_ext line_profiler
%autoreload
from   MassTodonPy  import MassTodon
import cPickle      as     pickle
from   MassTodonPy.MatchMaker   import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from   collections import Counter

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7


# modifications = {   ('N',2) :       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )
# path  = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
# path += 'Ubiquitin_ETD_10 ms_1071.mzXML'

file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.mzXML'
cutOff = 100; topPercent = .999

M.readSpectrum(path=file_path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)

M.prepare_problems(M_minProb)
mu=1e-5; lam=0.0; nu=0.001

# %%time
# MassTodonResults = M.run(solver='sequential', method='MSE', mu=mu, lam=lam, nu=0.001)
#
results_file = "/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.matteo"
# pickle.dump( MassTodonResults, open(results_file, "wb" ) )
with open(results_file, 'rb') as f:
    MassTodonResults = pickle.load(f)


%%time
Basic = reaction_analist_basic(MassTodonResults, Q, fasta)

%%time
Intermediate = reaction_analist_intermediate(MassTodonResults, Q, fasta)

%%time
UpperIntermediate = reaction_analist_upper_intermediate(MassTodonResults, Q, fasta, mu=0.0, lam=0.01)


print '\tETnoD\nBasic:', Basic['ETnoD'], '\nIntermediate:',Intermediate['ETnoD'], ' \nUpperIntermediate:',UpperIntermediate['ETnoD']
