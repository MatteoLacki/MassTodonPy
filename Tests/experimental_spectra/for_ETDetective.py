import cPickle as pickle
from    collections import Counter, defaultdict
import  pandas as pd


results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS/'

with open( results_path+'all_real_results', 'r' ) as handle:
    results = pickle.load(handle)

def get_name(key):
    return "_".join(map(str,key))

def get_subsequence(fasta, name):
    if name[0]=='p':
        return '*'+fasta+'*'
    if name[0]=='z':
        return fasta[len(fasta) - int(name[1:]):]+'*'
    else:
        return '*'+fasta[0:int(name[1:])]

def gen_Ciachable_data(R, fasta, Q):
    for r in R:
        name = r['molType']
        if name == 'precursor' and r['q'] == Q:
            f = {'seq':get_subsequence(fasta, name), 'Q':r['q'],'G':r['g'],'fragName':name, 'intensity':r['estimate'] }
            yield f
    for r in R:
        name = r['molType']
        if not (name == 'precursor' and r['q'] == Q):
            f = {'seq':get_subsequence(fasta, name), 'Q':r['q'],'G':r['g'],'fragName':name, 'intensity':r['estimate'] }
            yield f

directory = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/results/'


for params, Res, WH, WV, RA in results:
    fasta, Q, WH, WV, L, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha = params
    file_name = get_name(('seq-'+fasta,'WH-'+str(WH),'WV-'+str(WV)))
    R = []
    for res in Res:
        for a in res['alphas']:
            if a['estimate'] > 0:
                R.append(a)
    CiachableData = pd.DataFrame(gen_Ciachable_data(R, fasta, Q))[[4,1,0,2,3]]
    CiachableData.to_csv( path_or_buf = directory+file_name+'.csv', index = False, sep='\t')
