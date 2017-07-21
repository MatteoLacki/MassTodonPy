import  os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy import MassTodonize
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
import cPickle as pickle

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)

results = []

ubiquitin = ubiquitins[0]
for ubiquitin in ubiquitins:
    fasta= ubiquitin['fasta']
    Q    = ubiquitin['Q']
    res  = MassTodonize(fasta           = fasta,
                        precursor_charge= Q,
                        mz_prec         = .05,
                        joint_probability_of_envelope = .999,
                        spectrum        = ubiquitin['spectrum'],
                        opt_P           = .95,
                        solver          = 'multiprocessing',
                        multiprocesses_No = None,
                        max_times_solve = 10,
                        raw_data        = True,
                        verbose         = True )
    ubiquitin['experimental_setting']['Q'] = Q
    results.append({
        'masstodon_output':     res,
        'experimental_setting': ubiquitin['experimental_setting'],
        'etdetective_input':    results_to_etdetective( res, fasta, Q )
    })

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_only_real/ubiquitins.masstodon', 'w') as h:
    pickle.dump(results, h, protocol=pickle.HIGHEST_PROTOCOL)
