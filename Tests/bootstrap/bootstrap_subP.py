import  os
os.environ['OMP_NUM_THREADS'] = "1"

from    MassTodonPy.TestScripts.substancesP import substancesP
from    bootstrap_misc  import analyze_experiments
import  sys

_, results_path = sys.argv

if not os.path.exists(results_path):
    os.makedirs(results_path)


for mol in substancesP:
    mol['Q'] = mol['precursorCharge']
    del mol['precursorCharge']

# analyze_experiments(substancesP,
#                     results_path,
#                     K=2,
#                     mz_prec = .065,
#                     opt_P   = .99,
#                     bootstrap_size=10,
#                     verbose = True )

analyze_experiments(substancesP,
                    results_path,
                    K=None,
                    mz_prec = .065,
                    opt_P   = .99,
                    bootstrap_size=250,
                    verbose = False )
