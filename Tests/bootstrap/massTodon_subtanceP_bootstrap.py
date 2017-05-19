from run_massTodon_on_WH_WV import getResults
from collections import Counter
import numpy as np
import numpy.random as npr
import cPickle as pickle
from    read_experiments import experiments as substancesP

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

ionsNo, repetsNo = 100000, 2

for i,exp in enumerate(substancesP):
    fasta, Q, WH, WV, L, modifications, (Ms, Is_real) = exp
    R = MassTodon_bootstrap(exp, ionsNo, repetsNo)
    with open( results_path+'WH-'+str(WH)+'_WV-'+str(WV)+'_ID-'+str(i)+'.sadoMasto', 'w' ) as handle:
        pickle.dump(R, handle)
    print
    print 'Dumped',i,'out of',len(substancesP),'.'
    print
