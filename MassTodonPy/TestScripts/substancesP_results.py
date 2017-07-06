import pkg_resources
import cPickle as pickle

path = pkg_resources.resource_filename('MassTodonPy', 'Data/')

with open( path+'substancesP_results.example','r') as f:
    substancesP_results_macOS = pickle.load(f)
