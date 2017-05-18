from run_massTodon_on_WH_WV import getResults, data, parse_experiment
from collections import Counter, defaultdict
import numpy as np
import numpy.random as npr
import cPickle as pickle

experiments = [ parse_experiment(exp) for exp in data ]
# results     = [ getResults(*exp) for exp in experiments ]

results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/'
# with open( results_path+'substanceP_analyzed.sadoMasto', 'w' ) as handle:
#     pickle.dump(results, handle)

with open( results_path+'substanceP_analyzed.sadoMasto', 'r' ) as handle:
    results = pickle.load(handle)

R = defaultdict(list)
for params, Results, WH, WV, RA in results:
    R[(WH,WV)].append(RA)

R.keys()

R[(150,300)]
R[(125,300)]
