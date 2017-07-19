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

ubiquitin = ubiquitins[0]

ubis9 = [u for u in ubiquitins if u['precursorCharge'] == 9  and u['experimental_setting']['preActive'] == 'none' and u['experimental_setting']['supActive'] == 'OFF' and u['experimental_setting']['retentionTime'] == 2 and not u['experimental_setting']['files'].split('.')[0] in wrong_files ][0]
u['experimental_setting']['files']



from MassTodonPy import MassTodonize



res = MassTodonize( fasta           = ubis9['fasta'],
                    precursor_charge= ubis9['precursorCharge'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    spectrum        = ubis9['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    verbose         = False )

def iter():
    for l in res['raw_estimates']:
        for e in l['alphas']:
            if e['molType'] == 'precursor':
                yield e

list(iter())

res.keys()
2212994382 / float(2212994382 + 32353800)

res['basic_analysis']

Counter([ u['precursorCharge'] for u in ubiquitins if not u['experimental_setting'][u'files'].split('.')[0] in wrong_files])


###############################################################




import cPickle as pickle
import os
from collections import Counter
from pandas import DataFrame as DF
from itertools import islice

def make_a_row(info, real_or_bootstrap = 'real'):
    res = {'real_or_bootstrap':real_or_bootstrap}
    for k in ('ID','WV','WH'):
        res[k] = info[k]
    for algo in ('basic_analysis','intermediate_analysis','advanced_analysis'):
        for thing, thing_name in zip(info[algo], ('prob','count')):
            for k in thing:
                short_algo = algo.split('_')[0]
                res[short_algo+'-'+thing_name+'-'+str(k)] = thing[k]
    for k in info['summary']:
        res[k] = info['summary'][k]
    return res

def results_iter(res):
    yield make_a_row(res['real'], real_or_bootstrap = 'real')
    for info in res['bootstrap']:
        yield make_a_row(info, real_or_bootstrap = 'boot')



# indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_CSV_03_07_2017_mzPrec-065/'
indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_14_07_2017/'
outdir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_14_07_2017_csv/'


# outdir='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/CSV_03_07_2017_mzPrec-065/'

# if not os.path.exists(outdir):
    # os.makedirs(outdir)

need_to_break = False

def get_good():
    for r, d, fs in os.walk(indir):
        for f in islice(fs,None):
            print r+f
            with open(r+f,'r') as h:
                res = pickle.load(h)
            if not res['real']['files'].split('.')[0] in wrong_files:
                yield res
        # DF(results_iter(res)).to_csv(path_or_buf=outdir+f+str('.csv'), index=False)
res = next(get_good())
int(round(8500.0 / res['real']['precursorMZ']))



res['real']['retentionTime']
res['real']['preActive']
res['real']['supActive']
res['real'].keys()

res['real']['basic_analysis']
