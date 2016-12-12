%load_ext autoreload
%autoreload
from IsotopeCalculator import isotopeCalculator, atomCnt2string
import numpy as np
from collections import Counter, defaultdict
from math import exp, fsum


fasta ='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9
modifications = {}
ionsNo  = 1000000
P       = .999
isoCalc = isotopeCalculator()
digits  = 2
molecule= {'H':1000,'C':500,'O':200,'S':20,'N':200}

def agggregate( keys, values, digits=2):
    lists = defaultdict(list)
    keys  = keys.round(digits)
    for k, v in zip(keys, values):
        lists[k].append(v)
    newKeys     = np.array(lists.keys())
    newValues   = np.empty(len(newKeys))
    for k, v in zip( newKeys, np.nditer( newValues, op_flags=['readwrite'] )):
        v[...] = fsum(lists[k])
    return newKeys, newValues


# masses, probs = isoCalc.getEnvelope( molecule, P, digits)

masses, counts = isoCalc.randomSpectrum( molecule, Q, ionsNo, digits, P)
