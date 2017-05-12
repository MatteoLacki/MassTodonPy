import cPickle as pickle
from collections import Counter
import numpy as np
from math import sqrt
import pandas as pd
from MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate,  reaction_analist_upper_intermediate

fp_out = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/'
molsNo = 100000

with open(fp_out+'results_molsNo-'+str(molsNo), "rb") as f:
    res = pickle.load(f)

def get_correlation_and_L1_error(R):
    real, estimated = [], []
    for k in R:
        r,e = R[k]
        real.append(r)
        estimated.append(e)
    real, estimated = map(np.array,(real, estimated))
    corr = np.corrcoef(real, estimated)[0][1]
    L1error = sum(abs(R[k][0]-R[k][1]) for k in R) /( real.sum()+estimated.sum() )
    return corr, L1error

def analyze(res):
    for ((Q, fasta, eps, molsNo, probs), D),(R, FE, RFE, ra_res, res_tmp, simulation_res, spectrum) in res:
        corr, L1_error = get_correlation_and_L1_error(R)
        PTR, ETnoD, ETD= probs
        yield {'Q':Q,'eps':eps,'molsNo':molsNo,'PTR':PTR,'ETnoD':ETnoD,'ETD':ETD,'FE':FE,'RFE':RFE, 'corr':corr, 'L1':L1_error}

analyzed = pd.DataFrame(analyze(res))
analyzed.to_csv(path_or_buf=fp_out+'analyzed_'+str(molsNo)+'.csv', index=False)

# recalculating the parameters estimation ... correctly

def re_estimate_probs(res):
    i = 0
    Res = []
    for ((Q, fasta, eps, molsNo, probs), D),(R, FE, RFE, ra_res, res_tmp, simulation_res, spectrum) in res:
        di, mo = divmod(i,20)
        if mo==0:
            print i
        ra_res = {}
        PTR, ETnoD, ETD= probs
        try:
            ra_res['up_int']= reaction_analist_upper_intermediate(res_tmp, Q, fasta, 0.0) #last arg: minimal_estimated_intensity
        except:
            print 'Failed upper intermediate'
        try:
            ra_res['int']   = reaction_analist_intermediate(res_tmp, Q, fasta, 0.0)
        except:
            print 'Failed intermediate'
        try:
            ra_res['base']  = reaction_analist_basic(res_tmp, Q, fasta, 0.0)
        except:
            print 'Failed basic'
        i += 1
        Res.append(((Q,eps,molsNo,PTR,ETnoD,ETD,FE,RFE), ra_res))
    return Res

probs_reestimated = re_estimate_probs(res)

with open(fp_out+'results_molsNo-re-'+str(molsNo), "w") as f:
    pickle.dump(probs_reestimated,f)
