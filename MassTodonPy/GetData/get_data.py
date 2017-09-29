import pkg_resources
import cPickle as pickle

def get_data(dataset):
    assert dataset in ['substanceP', 'ubiquitin']
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open( path+dataset+'.example','r') as f:
        data = pickle.load(f)
    return data
