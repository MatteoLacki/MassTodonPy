import  json
import  cPickle as pickle

def change_key(seq, q, p, name):
    if name[0]=='c':
        name = 'c'+str(int(name[1:])-1)
    return name, q, p


# storagePath = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra.json'
# with open(storagePath) as data_file:
#     data = json.load(data_file)

with open('data/sigmas_probs.json', 'r') as f:
    s2p = json.load(f)

sigmas2probs = dict(s2p)
probs2sigmas = dict( (b,a) for a,b in s2p )
