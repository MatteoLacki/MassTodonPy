from run_massTodon_on_WH_WV import getResults, data, parse_experiment
from collections import Counter
import numpy as np
import numpy.random as npr

experiments = [ parse_experiment(exp) for exp in data ]
results = [ getResults(verbose=False, *exp) for exp in experiments ]

Counter(len(r) for r in results)
params, Results, WH, WV, RA = results[0]



def MassTodon_bootstrap(experiment, ionsNo, repetsNo, verbose=False):
    '''Perform bootstrap analysis of MassTodon results.'''
    fasta, Q, WH, WV, L, modifications, (Ms, Is_real) = experiment
    analisi = []
    for i, Is_sim in enumerate(npr.multinomial(ionsNo, Is_real/Is_real.sum(), repetsNo)):
        spectrum = (Ms, Is_sim)
        res = getResults(fasta, Q, WH, WV, L, modifications, spectrum, verbose=verbose)
        if len(res)>2:
            analisi.append(res)
        if divmod(i, 10)[1] == 0:
            print 'Finished',i+1,'out of',repetsNo,'.'
    return (fasta, Q, WH, WV, ionsNo, repetsNo), analisi

ionsNo, repetsNo = 100000, 100
R = MassTodon_bootstrap(experiments[0], ionsNo, repetsNo)
