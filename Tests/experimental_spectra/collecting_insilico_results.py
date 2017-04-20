import cPickle as pickle
import os
import pandas as pd
from collections import Counter
directory = '/Users/matteo/Documents/MassTodon/in_silico_results/masstodon_insilico/'

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
                totalE = sum(res[r] for r in res)
                totalR = sum(res[r] for r in mols)
                E = []
                R = []
                if totalR>0 and totalE>0:
                    for k in set(res) | set(mols):
                        R.append(mols[k]/float(totalR))
                        E.append(res[k]/float(totalE))
                    yield (sigma, molsNo, totalR), (R, E)

non_strange_results = list(iter_res(directory))
import json

# non_strange_results[100]

path='/Users/matteo/Documents/MassTodon/in_silico_results/wloczykij_latest.json'
with open(path, 'w') as fp:
    json.dump(non_strange_results, fp)
