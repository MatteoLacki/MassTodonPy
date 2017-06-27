import  cPickle         as pickle
from    bootstrap_misc  import analyze_experiments
import  sys
import  os

_, data_path, results_path = sys.argv
if not os.path.exists(data_path):
    print 'Fuck'

if not os.path.exists(results_path):
    os.makedirs(results_path)
with open(data_path+'substanceP_spectra_parsed.cPickle', 'r') as f:
    substancesP = pickle.load(f)

analyze_experiments(substancesP, results_path, K=2, bootstrap_size=10)
