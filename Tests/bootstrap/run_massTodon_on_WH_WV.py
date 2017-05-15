import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
from    MassTodonPy.MatchMaker   import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate

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



# fasta, Q, WH, WV, L, modifications, spectrum = exp
# cutOff = 100; topPercent = .999; max_times_solve=30
# jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# L1_x = L2_x = L1_alpha = L2_alpha = .001
# verbose = True
def getResults(fasta, Q, WH, WV, L, modifications, spectrum, jP=.999, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 100, topPercent = .999, max_times_solve=30, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False):
    params = (fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha)
    try:
        M = MassTodon(  fasta           = fasta,
                        precursorCharge = Q,
                        precDigits      = precDigits,
                        jointProbability= jP,
                        mzPrec          = mzPrec  )

        M.readSpectrum( spectrum        = spectrum,
                        cutOff          = cutOff,
                        digits          = precDigits,
                        topPercent      = topPercent  )

        M.prepare_problems(M_minProb)
        T0_deconv = time()
        Results = M.run(solver  = 'sequential',
                        method  = 'MSE',
                        max_times_solve = max_times_solve,
                        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                        verbose = verbose )
        T1_deconv = time()
        T_deconv  = T1_deconv - T0_deconv
        RA = {}

        try:
            RA['base'] = reaction_analist_basic(Results, Q, fasta, 0.0)
        except:
            print 'Base missing'
        try:
            RA['inter'] = reaction_analist_intermediate(Results, Q, fasta, 0.0)
        except:
            print 'Intermediate missing'
        try:
            RA['up_inter']  = reaction_analist_upper_intermediate(Results, Q, fasta, 0.0)
        except:
            print 'Upper Intermediate missing'

        res = (params, Results, WH, WV, RA)
    except Exception as e:
        res = (params, e)
    return res
