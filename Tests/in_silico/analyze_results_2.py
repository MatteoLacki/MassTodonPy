import cPickle as pickle
from pandas import DataFrame as DF
import os
from itertools import islice

indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/'
N = 1
outdir='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/wynyky.matteo'

def make_a_row(res):
    r = {'mols_no': res['molsNo'], 'sigma': res['sigma'] }
    r.update(res['summary'])
    r.update(res['probs'])
    for algo in ('basic_analysis','intermediate_analysis','advanced_analysis'):
        for thing, thing_name in zip(res[algo], ('prob','count')):
            for k in thing:
                short_algo = algo.split('_')[0]
                r[short_algo+'-'+thing_name+'-'+str(k)] = thing[k]
    for k in res['deconvolution_stats']:
        r['deconv_'+k] = res['deconvolution_stats'][k]
    return r


def results_iter(indir, N=None):
    for r, d, fs in os.walk(indir):
        for f in islice(fs,N):
            with open(r+f,'r') as h:
                res = pickle.load(h)
            yield make_a_row(res)

DF(results_iter(indir)).to_csv(path_or_buf=outdir, index=False)

res['sigma']
res['summary']
