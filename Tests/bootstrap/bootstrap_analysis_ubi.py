import cPickle as pickle
import os
from collections import Counter
from pandas import DataFrame as DF
from itertools import islice


def make_a_row(info, real_or_bootstrap = 'real'):
    res = {'real_or_bootstrap':real_or_bootstrap}
    for k in ('preActive', 'supActive', 'retentionTime', 'precursorMZ', 'ID', 'run', 'files'):
        if k in info:
            res[k] = info[k]
    for algo in ('basic_analysis','intermediate_analysis','advanced_analysis'):
        for thing, thing_name in zip(info[algo], ('prob','count')):
            for k in thing:
                short_algo = algo.split('_')[0]
                res[short_algo+'-'+thing_name+'-'+str(k)] = thing[k]
    for k in info['summary']:
        res[k] = info['summary'][k]
    return res

NoneTypesCnt = 0
def results_iter(res):
    if res['real'] != None:
        yield make_a_row(res['real'], real_or_bootstrap = 'real')
    for info in res['bootstrap']:
        if info != None:
            yield make_a_row(info, real_or_bootstrap = 'boot')



# indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/UBI/'
# outdir='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/UBI_csv/'

    # The most thresh data.
indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_14_07_2017/'
outdir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_14_07_2017_csv/'

# indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Boot_ubi_test/'
# outdir='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Boot_ubi_test_res/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

N = None
for r, d, fs in os.walk(indir):
    for f in islice(fs,N):
        with open(r+f,'r') as h:
            res = pickle.load(h)
        try:
            DF(results_iter(res)).to_csv(path_or_buf=outdir+f+str('.csv'), index=False)
        except:
            print f
            print 'Something went terribly wrong. Please contact the God, as he is ultimately responsible for it. '
            break
