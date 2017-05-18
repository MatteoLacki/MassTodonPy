import  json
import  numpy as np
import  cPickle as pickle
from    MassTodonPy  import MassTodon
from    time import time

storagePath = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/spectra.json'
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
    info    = (fasta, Q, WH, WV, L, modifications, spectrum)
    return info

experiments = [ parse_experiment(exp) for exp in data ]

def get_test_deconvolution_results(specNo = 0,
        jP=.999, mzPrec=.05, precDigits=2, M_minProb=.7,
        cutOff = 100., cutOff2=0.0, topPercent = .999, max_times_solve=30,
        L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False
    ):
    substanceP = not specNo is None
    if substanceP:
        fasta, Q, WH, WV, L, modifications, spectrum = experiments[specNo]
    else:
        fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
        Q=8
        file_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/FRL_220715_ubi_952_ETD_40ms_04.mzXML'
    M = MassTodon(  fasta           = fasta,
                    precursorCharge = Q,
                    precDigits      = precDigits,
                    jointProbability= jP,
                    mzPrec          = mzPrec  )

    if substanceP:
        M.readSpectrum( spectrum        = spectrum,
                        cutOff          = cutOff,
                        digits          = precDigits,
                        topPercent      = topPercent    )
    else:
        M.readSpectrum( path        = file_path,
                        cutOff      = cutOff,
                        digits      = precDigits,
                        topPercent  = topPercent    )
    M.prepare_problems(M_minProb)
    T0_deconv = time()
    Results = M.run(solver  = 'sequential',
                    method  = 'MSE',
                    max_times_solve = max_times_solve,
                    L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                    verbose = verbose )
    T1_deconv = time()
    T_deconv  = T1_deconv - T0_deconv
    return Results, T_deconv

subPres = get_test_deconvolution_results(specNo = 0, verbose=True)
ubiRes  = get_test_deconvolution_results(specNo = None, verbose=True)

results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/'
with open( results_path+'substanceP.sadoMasto', 'w' ) as handle:
    pickle.dump(subPres, handle)

with open( results_path+'ubiquitin.sadoMasto', 'w' ) as handle:
    pickle.dump(ubiRes, handle)
