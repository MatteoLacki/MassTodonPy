from MassTodonPy import MassTodonize
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from MassTodonPy.TestScripts.substanceP import substanceP
from time  import time
from pprint import pprint
import json

mol = substanceP.copy()
# mol = ubiquitin.copy()

# import cPickle as pickle
# from pandas import DataFrame as DF
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Data/ubiquitins.example', 'r') as h:
#     ubiquitins = pickle.load(h)

# mol = ubiquitins[120]

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    highcharts      = True,
                    verbose         = False )

highcharts = res['highcharts']



with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/highcharts.json', 'w') as f:
    json.dump(highcharts, f)




# import json
# H = make_etnod_ptr_probability_plot(algos)


# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_data.json','w') as h:
#     json.dump(H.data, h)
