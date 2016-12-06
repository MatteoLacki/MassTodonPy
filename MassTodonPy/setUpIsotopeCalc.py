%load_ext autoreload
%autoreload

# from Formulator import makeFragments
# from itertools  import chain
from IsotopeCalculator import isotopeCalculator
# , atomCnt2string
import numpy as np


fasta ='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9
modifications = {}
ionsNo  = 1000000
P       = .999

# prec, c, z = makeFragments(fasta, type='cz', modifications={})
# prec, c, z = [ list(s()) for s in [prec,c,z] ]

# atomCnt = prec[0]['atomCnt']
# atomCnt2string(atomCnt)

isoCalc = isotopeCalculator()
# masses, probs = isoCalc.isoEnvelope({'H':10,'O':2},P)
# isoCalc.isotopicEnvelopes.keys()
#
# %%time
# masses, probs = isoCalc.isoEnvelope({'H':1000,'O':200,'S':20},P)

# %%time
# masses, probs = isoCalc.isoEnvelope({'H':1000,'O':200,'S':20},P)
# isoCalc.isotopicEnvelopes.keys()

isoCalc.randomSpectrum()

from Formulator import makeFragments
