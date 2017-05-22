import cPickle as pickle
import os
import pandas as pd
from collections import Counter


analysis_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/'
results_path = analysis_path + 'RESULTS/'
ionsNo, repetsNo = 100000, 250
bootstrap = []

# file_name = 'WH-100_WV-300_ID-11.sadoMasto'
for file_name in os.listdir(results_path):
    with open( results_path+file_name, 'r' ) as f:
        R = pickle.load(f)
    if file_name.endswith(".sadoMasto"):
        ID = int(file_name.split('_')[2].split('.')[0][3:])
        R = R[0], R[1], ID
        bootstrap.append(R)
    else:
        real_data = R

def prepare_row(RA, WH, WV, algo, real_or_sim, count_or_prob, ID):
    row = {}
    row['algo'] = algo
    row['real_or_sim'] = real_or_sim
    row['WH'] = WH
    row['WV'] = WV
    row['ID'] = ID
    row['count_or_prob'] = count_or_prob
    if count_or_prob == 'prob':
        idx_cop = 0
    else:
        idx_cop = 1
    row.update( dict(RA[algo][idx_cop]) ) # this will give probs
    return row

def get_real(real_data):
    ID = 0
    for params, Results, WH, WV, RA in real_data:
        for algo in RA:
            yield prepare_row(RA, WH, WV, algo, 'real', 'prob', ID)
            yield prepare_row(RA, WH, WV, algo, 'real', 'count', ID)
        ID += 1

def get_sim(bootstrap):
    for (fasta, Q, WH, WV, molsNo, repetsNo), sims, ID in bootstrap:
        for RA in sims:
            for algo in RA:
                yield prepare_row(RA, WH, WV, algo, 'sim', 'prob', ID)
                yield prepare_row(RA, WH, WV, algo, 'sim', 'count', ID)

D = pd.DataFrame(get_real(real_data))
S = pd.DataFrame(get_sim(bootstrap))

D.to_csv( path_or_buf = analysis_path+'real_data_3.csv', index = False)
S.to_csv( path_or_buf = analysis_path+'simulations_3.csv', index = False)
