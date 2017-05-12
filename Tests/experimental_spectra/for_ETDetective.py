import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
import  cPickle      as     pickle
from    MassTodonPy.MatchMaker   import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from    collections import Counter, defaultdict
import  pandas as pd

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
# cutOff = 100; topPercent = .999
# jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# mu=1e-5; lam=0.0; nu=0.001
def getResults(fasta, Q, WH, WV, L, modifications, spectrum, jP=.999, mzPrec=.05, precDigits=2, M_minProb=.7, cutOff = 100, topPercent = .999, max_times_solve=30, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001):
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
                        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha )
        T1_deconv = time()
        T_deconv  = T1_deconv - T0_deconv

        T0_basic = time()
        Basic = reaction_analist_basic(Results, Q, fasta, 0.0)
        T1_basic = time()
        T_basic  = T1_basic-T0_basic

        T0_inter = time()
        Intermediate = reaction_analist_intermediate(Results, Q, fasta, 0.0)
        T1_inter = time()
        T_inter  = T1_inter-T0_inter

        T0_up_inter = time()
        UpperIntermediate = reaction_analist_upper_intermediate(Results, Q, fasta, 0.0)
        T1_up_inter = time()
        T_up_inter = T1_up_inter-T0_up_inter

        optimal, nonoptimal, totalError = M.flatten_results()

        res = (True, params, Results, T_deconv, Basic, T_basic, Intermediate, T_inter, UpperIntermediate, T_up_inter, optimal, nonoptimal, totalError, WH, WV)

    except Exception as e:
        res = (False, params, e)
    return res

experiments = [ parse_experiment(exp) for exp in data ]
results = [ getResults(*exp) for exp in experiments ]

def get_name(key):
    return "_".join(map(str,key))

def get_subsequence(fasta, name):
    if name[0]=='p':
        return '*'+fasta
    if name[0]=='z':
        return fasta[len(fasta) - int(name[1:]):]
    else:
        return fasta[0:len(fasta)]

def gen_Ciachable_data(R, fasta):
    for r in R:
        name = r['molType']
        f = {'seq':get_subsequence(fasta, name), 'Q':r['q'],'G':r['g'],'fragName':name, 'intensity':r['estimate'] }
        yield f

directory = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/results/'

for r in results:
    fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha = r[1]
    file_name = get_name(('seq-'+fasta,'WH-'+str(WH),'WV-'+str(WV)))
    R = []
    for res in r[2]:
        for a in res['alphas']:
            if a['estimate'] > 0:
                R.append(a)
    CiachableData = pd.DataFrame(gen_Ciachable_data(R, fasta))[[4,1,0,2,3]]
    CiachableData.to_csv( path_or_buf = directory+file_name+'.csv', index = False)
