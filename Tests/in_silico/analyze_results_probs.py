import cPickle as pickle
import pandas as pd
from collections import Counter, defaultdict
import numpy as np
from itertools import product


fp_out = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/'
molsNo = 1000

with open(fp_out+'results_molsNo-re-'+str(molsNo), "rb") as f:
    probs = pickle.load(f)

def get_prob(ra, ra_type, key):
    if ra_type=='base':
        return ra[ra_type].get(key,None)
    else:
        return ra[ra_type][0].get(key,None)


def gen_probs(probs):
    for (Q,eps,molsNo,PTR,ETnoD,ETD,FE,RFE), ra_res in probs:
        P = {}
        for ra_type, key in product( ('base','int','up_int'),('PTR','ETnoD','fragmentation','reaction') ):
            P[key+'_'+ra_type] = get_prob(ra_res, ra_type, key)
        info = {'Q':Q,'eps':eps,'molsNo':molsNo,'PTR':PTR,'ETnoD':ETnoD,'ETD':ETD,'FE':FE,'RFE':RFE}
        P.update(info)
        yield P

Res = pd.DataFrame(gen_probs(probs))
Res.to_csv(path_or_buf=fp_out+'probs_'+str(molsNo)+'.csv', index=False)
