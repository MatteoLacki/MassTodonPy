import json
import numpy as np
import cPickle as pickle

with open( '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubi_spectra.json', 'r') as h:
    spectra = json.load(h)

with open( '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubi_params.json', 'r') as h:
    params = json.load(h)


print 'The same data:', len(spectra) == len(params)

ubiquitin_fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

ubiquitins = []
for spectrum, param in zip(spectra, params):
    Q = param['precursorCharge']
    del param['precursorCharge']
    mol = {'fasta':    ubiquitin_fasta,
           'modifications':         {},
           'experimental_setting':  param,
           'precursorCharge':       Q }
    mzs         = np.empty(len(spectrum))
    intensities = np.empty(len(spectrum))
    for i, dic in enumerate(spectrum):
        mzs[i] = dic['mz']
        intensities[i] = dic['intensity']
    mol['spectrum'] = mzs, intensities
    ubiquitins.append(mol)

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/ubiquitins.example', 'w') as h:
    pickle.dump(ubiquitins, h)
