import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
import  cPickle      as     pickle
from    collections import Counter, defaultdict

storagePath = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/substanceP_spectra.json'

with open(storagePath) as data_file:
    data = json.load(data_file)

def extract_WH_WV(s):
    s = str(s)
    if ' ' in s:
        s = s.replace('1,5','150').replace(' ','')
    WH, WV = s.split('-')[-2:]
    WH = int(WH[2:])
    WV = int(WV[2:])
    return WH, WV

def parseSpectrum(exp):
    N = len(exp['mass_spectrum'])
    mz = np.empty(N)
    intensities = np.empty(N)
    for n in xrange(N):
        mz[n] = exp['mass_spectrum'][n]['mz']
        intensities[n] = exp['mass_spectrum'][n]['intensity']
    return mz, intensities

def parseTerminus(terminus):
    return dict([ (str(atom), terminus[atom][0]) for atom in terminus])

def parse_experiment(exp):
    fasta   = str(exp['precursor'][0])
    Q       = exp['maxProtonsNo'][0]
    WH, WV  = extract_WH_WV(exp['instrumental_setting'][0])
    L = len(fasta)
    modifications = {}
    for terminus in set(exp.keys()) & set(['cTerminus', 'nTerminus']):
        if terminus[0] == 'c':
            key = ('N',1)
        else:
            key = ('C',L)
            modifications[key] = parseTerminus(exp[terminus])
    spectrum= parseSpectrum(exp)
    info    = fasta, Q, WH, WV, L, modifications, spectrum
    return info

experiments = [ parse_experiment(exp) for exp in data ]
fasta, Q, WH, WV, L, modifications, spectrum = experiments[0]
modifications = { 'C11': {'H':1,'O':-1,'N':1} }

substanceP = {  'name'      : 'substanceP',
                'fasta'     : fasta,
                'Q'         : Q,
                'WH'        : WH,
                'WV'        : WV,
                'modifications': modifications,
                'spectrum'  : spectrum  }

standard_file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/substanceP.example'

with open(standard_file_path,'w') as f:
    pickle.dump(substanceP, f)
