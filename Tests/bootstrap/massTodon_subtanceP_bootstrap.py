from run_massTodon_on_WH_WV import getResults, data, parse_experiment
from collections import Counter
import numpy as np
import numpy.random as npr
import cPickle as pickle

experiments = [ parse_experiment(exp) for exp in data ]
results     = [ getResults(*exp) for exp in experiments ]


def MassTodon_bootstrap(experiment, ionsNo, repetsNo, verbose=False):
    '''Perform bootstrap analysis of MassTodon results.'''
    fasta, Q, WH, WV, L, modifications, (Ms, Is_real) = experiment
    analisi = []
    for i, Is_sim in enumerate(npr.multinomial(ionsNo, Is_real/Is_real.sum(), repetsNo)):
        spectrum = (Ms, Is_sim)
        res = getResults(fasta, Q, WH, WV, L, modifications, spectrum, verbose=verbose)
        if len(res)==5:
            analisi.append(res[4])
        if divmod(i, 10)[1] == 0:
            print 'Finished',i+1,'out of',repetsNo,'.'
    return (fasta, Q, WH, WV, ionsNo, repetsNo), analisi

results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS/'

with open( results_path+'all_real_results', 'w' ) as handle:
    pickle.dump(results, handle)

ionsNo, repetsNo = 100000, 250

for i,exp in enumerate(experiments):
    fasta, Q, WH, WV, L, modifications, (Ms, Is_real) = exp
    R = MassTodon_bootstrap(exp, ionsNo, repetsNo)
    with open( results_path+'WH-'+str(WH)+'_WV-'+str(WV)+'_ID-'+str(i)+'.sadoMasto', 'w' ) as handle:
        pickle.dump(R, handle)
    print
    print 'Dumped',i,'out of',len(experiments),'.'
    print
