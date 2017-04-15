import cPickle as pickle

with open('test.matteo', 'w') as f:
    pickle.dump( 'A', f)

with open('test.matteo', 'w') as f:
    pickle.dump( 'B', f)

with open('test.matteo', 'rb') as f:
    x = pickle.load(f)

x
