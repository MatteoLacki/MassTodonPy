import cPickle as pickle

path = '/Users/matteo/Wloczykij/masstodon/masstodon_insilico/sigma_999999426697_molsNo_1.matteo'

with open(path, 'rb') as f:
    x = pickle.load(f)

x
