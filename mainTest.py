%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules, makeFragments, pandizeSubstances
from Formulator.fasta2atomcnt import fasta2atomCnt
from InSilico.spectrumGenerator import insilicoSpectrum, genIsotopicEnvelope

ubiquitin   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
substanceP  = 'RPKPQQFFGLM'
fastas      = [substanceP, ubiquitin, ubiquitin+ubiquitin+ubiquitin]
fasta = substanceP
Q = 3
# modifications = {   ('N',2) :       {'H': -1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': -1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': -2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': -2, 'S': +2, 'N': +2} }
modifications = {}
ionsNo = 100000
P = .999

# data = pandizeSubstances(*makeFragments(fasta, 'cz', modifications))
# molecules = list( genMolecules(fasta, 3, 'cz', modifications) )
# fasta2atomCnt(fasta)

IS = insilicoSpectrum(fasta, Q, ionsNo, P)
sample = IS.rvs(ionsNo)
len(sample)
