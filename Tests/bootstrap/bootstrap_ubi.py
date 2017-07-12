import  os
os.environ['OMP_NUM_THREADS'] = "1"

import  sys
from    bootstrap_misc  import analyze_experiments
import cPickle as pickle
with open('/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)


_, results_path = sys.argv

if not os.path.exists(results_path):
    os.makedirs(results_path)

analyze_experiments(ubiquitins,
                    results_path,
                    # K = 2,
                    mz_prec         = .05,
                    # cut_off         = 500, # this is totally nonsensical
                    opt_P           = .95, # cannot take .8 / cannot take .99
                    bootstrap_size  = 250,
                    verbose         = True     )
