from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
import json
from pandas import DataFrame as DF

import cPickle as pickle
from pandas import DataFrame as DF

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)


# mol = substanceP.copy()
mol = ubiquitins[120]

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['precursorCharge'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    raw_data = True,
                    for_plot = True,
                    verbose  = False )


with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_original.json', 'w') as handler:
    json.dump(zip(*res['spectra']['original']), handler)


long_data = res['for_plot_long']
DF(long_data).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data.csv',
    index = False )

# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data.json', 'w') as handler:
#     json.dump(long_data, handler)
