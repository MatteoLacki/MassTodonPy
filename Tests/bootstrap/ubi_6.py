import  os
os.environ['OMP_NUM_THREADS'] = "1"

import  sys
from    bootstrap_misc  import analyze_experiments
import  cPickle as pickle
from collections import Counter

###############################################################
initial_run_path = '/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun'
wrong_files = [f.split('.')[0] for f in list(os.walk(initial_run_path))[0][2]]
###############################################################


# with open('../../MassTodonPy/Data/ubiquitins.example', 'r') as h:
#     ubiquitins = pickle.load(h)

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)



ubis6 = [u for u in ubiquitins if u['precursorCharge'] == 6  and u['experimental_setting']['preActive'] == 'none' and u['experimental_setting']['supActive'] == 'OFF' and u['experimental_setting']['retentionTime'] == 2 and not u['experimental_setting']['files'].split('.')[0] in wrong_files ][0]

ubis6['experimental_setting']['files']

from MassTodonPy import MassTodonize


res05 = MassTodonize( fasta         = ubis6['fasta'],
                    precursor_charge= ubis6['precursorCharge'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    spectrum        = ubis6['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    verbose         = False )

res01 = MassTodonize( fasta         = ubis6['fasta'],
                    precursor_charge= ubis6['precursorCharge'],
                    mz_prec         = .01,
                    joint_probability_of_envelope = .999,
                    spectrum        = ubis6['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    verbose         = False )

# mass = ubis6['spectrum'][0]
# intensity = ubis6['spectrum'][1]
from pandas import DataFrame as DF

# DF(mass, intensity).to_csv(path_or_buf='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/dupa.csv', index=True)

def iter(res):
    for l in res['raw_estimates']:
        for e in l['alphas']:
            if e['molType'] == 'precursor':
                yield e

d05 = list(iter(res05))
d05


d01 = list(iter(res01))



# ETnoD_cnt = sum( e['estimate']*e['g'] for e in d )
# PTR_cnt = sum( e['estimate']*( Q - e['q'] - e['g']) for e in d )
# ETnoD_cnt/float(ETnoD_cnt + PTR_cnt)

res['basic_analysis'][0]['ETnoD']
