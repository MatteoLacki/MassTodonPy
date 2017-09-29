from MassTodonPy import MassTodonize

import pkg_resources
path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
import cPickle as pickle
with open( path+'substanceP.example','r') as f:
    substanceP = pickle.load(f)


from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time

mol = substanceP.copy()
# mol = ubiquitin.copy()

res = MassTodonize( fasta                           = mol['fasta'],
                    precursor_charge                = mol['Q'],
                    mz_prec                         = .05,
                    joint_probability_of_envelope   = .999,
                    modifications                   = mol['modifications'],
                    spectrum                        = mol['spectrum'],
                    opt_P                           = .99,
                    solver                          = 'multiprocessing',
                    multiprocesses_No               = None,
                    distance_charges                = 5,
                    max_times_solve                 = 10,
                    raw_data                        = True,
                    output_csv_path                 = '/Users/matteo/Documents/MassTodon/results/',
                    highcharts                      = False,
                    verbose                         = False )

res['summary']
