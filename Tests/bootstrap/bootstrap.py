import  os
os.environ['OMP_NUM_THREADS'] = "1"

from    MassTodonPy.TestScripts.standard_datasets import substancesP
from    bootstrap_misc  import analyze_experiments
import  sys

_, results_path = sys.argv

if not os.path.exists(results_path):
    os.makedirs(results_path)

# analyze_experiments(substancesP, results_path, K=2, bootstrap_size=200)
# analyze_experiments(
#     substances      = substancesP,
#     results_path    = results_path,
#     K = 2, bootstrap_size = 10, mz_prec = .065, verbose = True )
# analyze_experiments(substancesP, results_path, bootstrap_size=250)

analyze_experiments(    substancesP,
                        results_path,
                        mz_prec         = .065,
                        opt_P           = .99,
                        bootstrap_size  = 250   )
