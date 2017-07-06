import pkg_resources
import cPickle as pickle

path = pkg_resources.resource_filename('MassTodonPy', 'Data/')

with open( path+'substancesP.example','r') as f:
    substancesP = pickle.load(f)
