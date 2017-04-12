from IsoSpecPy import IsoSpecPy
from math import exp
import numpy as np
from math import exp
from collections import Counter

i = IsoSpecPy.IsoSpec.IsoFromFormula("H20000O10000", 0.9)
masses, logprobs, _ = i.getConfsRaw()

# for m,lp in zip(masses,logprobs):
#     print m, exp(lp)

massPrecDigits = 3
aggregator = Counter()
for mass, logprob in zip(masses, logprobs):
    aggregator[ round(mass, massPrecDigits) ] += exp(logprob)

masses = np.array(aggregator.keys())
probs  = np.empty(len(masses))
for i, mass in enumerate(masses):
    probs[i] = aggregator[mass]

q = 4
g = 5
masses
W = masses.copy()
W += 10
masses
W
masses = (masses + g)/(q+g)
masses.round(decimals=3)
