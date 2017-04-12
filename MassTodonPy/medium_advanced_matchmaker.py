%load_ext autoreload
%load_ext line_profiler
%autoreload
from   MassTodon    import MassTodon
import cPickle      as     pickle
from   MatchMaker   import reaction_analist_intermediate, reaction_analist_basic, reaction_analist_advanced
from   math         import exp

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

%%time
Intermediate = reaction_analist_intermediate(MassTodonResults, Q, fasta)

%%time
Basic = reaction_analist_basic(MassTodonResults, Q, fasta)

%%time
Advanced = reaction_analist_advanced(MassTodonResults, Q, fasta, maxIter=100, const=1000, eps = 0.0, crit='logLikDiff', verbose=False, tol=0.01)

Basic
Advanced
