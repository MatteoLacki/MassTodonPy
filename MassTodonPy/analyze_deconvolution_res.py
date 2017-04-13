import cPickle as pickle

path = '/Users/matteo/Documents/MassTodon/in_silico_results/wloczykij/sigma_3_molsNo_1.matteo'

with open(path, 'rb') as f:
    x = pickle.load(f)
