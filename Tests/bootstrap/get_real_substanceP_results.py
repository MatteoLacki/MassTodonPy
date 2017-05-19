from run_massTodon_on_WH_WV import getResults
from collections import Counter
import numpy as np
import numpy.random as npr
import cPickle as pickle
from    read_experiments import experiments as substancesP

results = [ getResults(*exp) for exp in substancesP ]
results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS/'
with open( results_path+'all_real_results', 'w' ) as handle:
    pickle.dump(results, handle)
