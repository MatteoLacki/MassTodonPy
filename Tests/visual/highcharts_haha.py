from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
from pprint import pprint
import json

# mol = substanceP.copy()
mol = ubiquitin.copy()

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


h1 = res['highcharts'][0]
h1.data
h1.options
h1.buildhtml()

for i, h in enumerate(res['highcharts']):
    h.save_file('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/h'+str(i))
    with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/d'+str(i)+'.json', 'w') as f:
        json.dump(h.data, f)



# import json
# H = make_etnod_ptr_probability_plot(algos)


# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_options.json','w') as h:
#     json.dump(H.option, h)
#
# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/etnodPtr_data.json','w') as h:
#     json.dump(H.data, h)
