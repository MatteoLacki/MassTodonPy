import  json
import  numpy as np
from    MassTodonPy  import MassTodon
import  cPickle      as     pickle
from    collections import Counter, defaultdict

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

cutOff = 100; topPercent = .999
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
mu=1e-5; lam=0.0; nu=0.001
Q = 8

M = MassTodon(  fasta           = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

path='/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/temp_data/'
path_file = path+'Ubiquitin_ETD_10 ms_1071.mzXML'

M.readSpectrum( path            = path_file,
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent  )

spectrum = M.spectrum

ubiquitin = {   'name'      : 'ubiquitin',
                'fasta'     : fasta,
                'Q'         : Q,
                'exp_param' : '10 ms',
                'modifications': {},
                'spectrum'  : spectrum  }

standard_file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/ubiquitin.example'

with open(standard_file_path,'w') as f:
    pickle.dump(ubiquitin, f)
