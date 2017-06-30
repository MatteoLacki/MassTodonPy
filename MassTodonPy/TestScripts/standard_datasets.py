import pkg_resources
import cPickle as pickle

path = pkg_resources.resource_filename('MassTodonPy', 'Data/')

with open( path+'ubiquitin.example','r') as f:
    ubiquitin = pickle.load(f)

with open( path+'substanceP.example','r') as f:
    substanceP = pickle.load(f)

with open( path+'substancesP.example','r') as f:
    substancesP = pickle.load(f)

with open( path+'substancesP_results.example','r') as f:
    substancesP_results_macOS = pickle.load(f)
