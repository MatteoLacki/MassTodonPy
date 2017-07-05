import cPickle as pickle
import os
from collections import Counter


# res = []
# for x in xrange(51):
#     with open('RESULTS_CSV_27_06_2017_czczmiel/'+str(x)) as f:
#         res.append(pickle.load(f))


# res = []
# for x in xrange(0):
#     with open('RESULTS_fixing_bootstrap/'+str(x)) as f:
#         res.append(pickle.load(f))

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_fixing_bootstrap/0','r') as f:
    res = pickle.load(f)

res = [res]

len(res[0]['real'])

len(res[0]['bootstrap'])


[b['summary'] for b in res[0]['bootstrap']]

errors = [ (i,len([ k for k in r['bootstrap'] if k is None])) for i,r in enumerate(res)  ]

i_max, error_max = max(errors, key=lambda x:x[1])

res[i_max].keys()


len(res[i_max]['bootstrap'])
res[i_max]['bootstrap'][1]
