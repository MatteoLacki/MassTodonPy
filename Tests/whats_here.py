import cPickle as pickle
import json

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/substanceP.sadoMasto', 'r' ) as f:
        R = pickle.load(f)

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra.json', 'r' ) as f:
        S = json.load(f)

len(S)

S[0].keys()

S[0]['instrumental_setting']
