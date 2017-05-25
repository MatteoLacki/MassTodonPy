import pkg_resources
import cPickle as pickle

path = pkg_resources.resource_filename('MassTodonPy', 'Data/')

with open( path+'ubiquitin.example','r') as f:
    ubiquitin = pickle.load(f)

with open( path+'substanceP.example','r') as f:
    substanceP = pickle.load(f)
