import cPickle as pickle
import os
import pandas as pd

directory = '/Users/matteo/Documents/MassTodon/in_silico_results/wloczykij/'

# filename = os.listdir(directory)[10]

def iter_res(directory):
    res = []
    for filename in os.listdir(directory):
        _,sigma,_,molsNo =  filename.replace('.matteo','').split('_')
        sigma = float('0.'+sigma)
        molsNo= int(molsNo)
        with open( directory + filename, 'rb') as f:
            res = pickle.load(f)
            for r in res:
                mols, x0, sigma, (error, non_optimal) = r
                yield {'molsNo':molsNo, 'x0':x0, 'sigma': sigma, 'error':error}

Res = pd.DataFrame(iter_res(directory))
Res.to_csv('/Users/matteo/Documents/MassTodon/in_silico_results/analysis.csv',
            index=False )
