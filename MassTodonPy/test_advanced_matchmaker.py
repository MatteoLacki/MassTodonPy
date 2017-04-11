%load_ext autoreload
%load_ext line_profiler
%autoreload
from   MassTodon  import MassTodon
from   MatchMaker import reaction_analist_advanced
from    math      import log, lgamma, log10, exp, sqrt
import cPickle as pickle

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

%%time
LogProb = reaction_analist_advanced(MassTodonResults, Q, fasta, maxIter=10, const=1000, eps = 0.0, verbose=False)

for r in LogProb:
    print r, exp(LogProb[r])
