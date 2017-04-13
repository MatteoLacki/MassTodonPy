%load_ext autoreload
%load_ext line_profiler
%autoreload
from   MassTodon    import MassTodon
import cPickle      as     pickle
from   MatchMaker   import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from    collections import Counter

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)


%%time
Basic = reaction_analist_basic(MassTodonResults, Q, fasta)

%%time
Intermediate = reaction_analist_intermediate(MassTodonResults, Q, fasta)

%%time
UpperIntermediate = reaction_analist_upper_intermediate(MassTodonResults, Q, fasta, mu=0.0, lam=0.01)

print '\tETnoD\nBasic:', Basic['ETnoD'], '\nIntermediate:',Intermediate['ETnoD'], ' \nUpperIntermediate:',UpperIntermediate['ETnoD']
