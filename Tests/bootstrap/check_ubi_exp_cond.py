import cPickle as pickle
import json
with open('/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)

len(ubiquitins)
ubiquitins[0]

mapping = {}
for ID, ubi in enumerate(ubiquitins):
    mapping[ID] = ubi['experimental_setting']['files']

indir = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/'
with open(indir+'ID_file.json', 'w') as h:
    json.dump(mapping, h)
