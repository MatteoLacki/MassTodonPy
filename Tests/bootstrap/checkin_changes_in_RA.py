from    MassTodonPy  import MassTodon
from    MassTodonPy.MatchMaker import czMatchMakerBasic, czMatchMakerIntermediate, czMatchMakerUpperIntermediate
import  cPickle as pickle
from    collections     import Counter

results_path= '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/'
substanceP  = False

if substanceP:
    with open( results_path+'substanceP.sadoMasto', 'r' ) as handle:
        Results, time = pickle.load(handle)
    Q=3
    fasta = 'RPKPQQFFGLM'
else:
    with open( results_path+'ubiquitin.sadoMasto', 'r' ) as handle:
        Results, time = pickle.load(handle)
    fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
    Q=8

basic = czMatchMakerBasic(Results, Q, fasta, accept_nonOptimalDeconv = False, min_acceptEstimIntensity = 100., verbose=False)
Probs, Counts = basic.pair()


intermediate = czMatchMakerIntermediate(Results, Q, fasta, accept_nonOptimalDeconv=False, min_acceptEstimIntensity=100., verbose=False)
Probs_inter, Counts_inter = intermediate.pair()

Counts_inter['ETnoD']
Counts_inter['ETnoD_frag']
Counts_inter['ETnoD_precursor']

Counts_inter['PTR_precursor']

Counts_inter['ETnoD_frag']
Counts_inter['ETnoD_precursor']


Counts_inter['PTR_frag']/float(Counts_inter['ETnoD_frag'])
Counts_inter['PTR_precursor']/float(Counts_inter['ETnoD_precursor'])

PTR = Counts_inter['PTR_precursor']+Counts_inter['PTR_frag']
ETnoD = Counts_inter['ETnoD_frag']+Counts_inter['ETnoD_precursor']

PTR/float(ETnoD)

upperIntermediate = czMatchMakerUpperIntermediate(Results, Q, fasta, accept_nonOptimalDeconv=False, min_acceptEstimIntensity=100., verbose=False)
upperIntermediate.pair()
