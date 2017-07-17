import cPickle as pickle
from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time import time
from math import log

mol = substanceP.copy()
# mol = ubiquitin.copy()

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .99,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    highcharts      = True,
                    verbose         = False )


results_4_etdetective = results_to_etdetective(res, mol['fasta'], mol['Q'], mol['modifications'])


with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/linking_to_ETDetective/test_input.json', 'w') as h:
    pickle.dump(results_4_etdetective, h)

# probs, counter = res['basic_analysis']
#
# -log(probs['no fragmentation'])/mol['Q']**2
