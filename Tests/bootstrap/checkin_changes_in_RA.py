from    MassTodonPy  import MassTodon
from    MassTodonPy.MatchMaker import czMatchMakerBasic, czMatchMakerIntermediate, czMatchMakerUpperIntermediate
import  cPickle as pickle

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
basic.pair()

intermediate = czMatchMakerIntermediate(Results, Q, fasta, accept_nonOptimalDeconv=False, min_acceptEstimIntensity=100., verbose=False)
intermediate.pair()

upperIntermediate = czMatchMakerUpperIntermediate(Results, Q, fasta, accept_nonOptimalDeconv=False, min_acceptEstimIntensity=100., verbose=False)
upperIntermediate.pair()
