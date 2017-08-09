import pkg_resources
import cPickle as pickle

if __name__ == '__main__':
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')

    with open( path+'substanceP.example','r') as f:
        substanceP = pickle.load(f)
