import  cPickle as pickle
import  os
import  pandas as pd
from    collections import Counter
import  json
directory = '/Users/matteo/Documents/MassTodon/in_silico_results/masstodon_insilico/'

probs = [   0.99999943,  0.99919849,  0.98831003,  0.95667342,  0.90814916,
            0.8518448 ,  0.79452072,  0.73981509,  0.68929589,  0.64340622,
            0.60204112,  0.56485823,  0.53143491,  0.50134371,  0.4741858 ,
            0.44960344,  0.4272823 ,  0.40694922,  0.388368  ,  0.37133488,
            0.35567401,  0.34123345,  0.32788163,  0.31550433,  0.30400205,
            0.29328792,  0.28328577,  0.27392868,  0.2651576 ,  0.25692031,
            0.24917047,  0.24186683,  0.2349726 ,  0.22845482,  0.22228395,
            0.21643338,  0.21087914,  0.20559956,  0.20057498,  0.19578757,
            0.19122111,  0.18686079,  0.18269311,  0.17870568,  0.17488715,
            0.17122712,  0.16771598,  0.16434491,  0.16110575,  0.15799095,
            0.15499353,  0.15210704,  0.14932545,  0.14664319,  0.14405505,
            0.14155621,  0.13914213,  0.1368086 ,  0.13455168,  0.13236767  ]


sigmas = set()
for filename in os.listdir(directory):
    _,sigma,_,molsNo = filename.replace('.matteo','').split('_')
    sigmas.add(float('0.'+sigma))

list(sigmas).sort()
Sigma2Prob = dict((sigma, 1.0 - prob) for sigma, prob in zip(sigmas, probs))

def iter_res(directory):
    res = []
    for filename in os.listdir(directory):
        _,sigma,_,molsNo = filename.replace('.matteo','').split('_')
        sigma = float('0.'+sigma)
        molsNo= int(molsNo)
        with open( directory + filename, 'rb') as f:
            results = pickle.load(f)
            # _, _, sigma, (res, mols) = results[-1]
            for _, _, sigma, (res, mols) in results:
                totalR = sum(mols[r] for r in mols)
                totalE = sum(res[r] for r in res)
                R = []
                E = []
                if totalR>0 and totalE>0:
                    for k in set(res) | set(mols):
                        R.append(mols[k]/float(totalR))
                        E.append(res[k]/float(totalE))
                    yield (sigma, Sigma2Prob[sigma],  molsNo, totalR), (R, E)

non_strange_results = list(iter_res(directory))
len(non_strange_results)
# non_strange_results[1768]

path = '/Users/matteo/Documents/MassTodon/in_silico_results/uniform.json'
with open(path, 'w') as fp:
    json.dump(non_strange_results, fp)
