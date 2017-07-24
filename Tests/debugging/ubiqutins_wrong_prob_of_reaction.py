from MassTodonPy import MassTodonize
from MassTodonPy.Formulator import make_formulas
import cPickle as pickle
from pandas import DataFrame as DF
from cvxopt import matrix
import numpy as np

# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubiquitins.example', 'r') as h:
#     ubiquitins = pickle.load(h)

file_name = u'Ubiquitin_ETD_952_25 ms'
with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/Ubiquitin_ETD_952_25_ms.example', 'r') as h:
    mol = pickle.load(h)


# file_name = u'Ubiquitin_ETD_1428_pre_act (NCE 20)_1 ms'
# The wrong spectrum.
# [ u['experimental_setting']['files'] for u in ubiquitins ]
# mol =[ u for u in ubiquitins if u['experimental_setting']['files'] == file_name + '.mzXML'][0]


# This is an initial run file!!!
res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    # precursor_charge= 8,
                    mz_prec         = .1,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    for_plot        = True,
                    verbose         = True )
                    # L1_x            = 1.0,
                    # L2_x            = 1.0,
                    # L1_alpha        = 1.0,
                    # L2_alpha        = 1.0 )


param = [ r for r in res['raw_estimates'] if r['status'] != 'optimal'][0]['param']

P = np.matrix(matrix(param['P']))
U, s, V = np.linalg.svd(P)
np.linalg.cond(P)


min(s)/max(s)
3.999840006355833e-05


param_ok = [ r for r in res['raw_estimates'] if r['status'] == 'optimal'][0]['param']
P_ok = np.matrix(matrix(param_ok['P']))
U, s, V = np.linalg.svd(P_ok)
min(s)*max(s)

np.linalg.cond(P_ok)


x = [ np.linalg.cond(np.matrix(matrix(r['param']['P']))) for r in res['raw_estimates'] if r['status'] == 'optimal']

[ r['sol']['iterations'] for r in res['raw_estimates'] if r['status'] != 'optimal' ]
y = [ r['sol']['iterations'] for r in res['raw_estimates'] if r['status'] == 'optimal' ]

import seaborn as sns
from pandas import DataFrame as DF


import matplotlib.pyplot as plt
# Show the results of a linear regression within each dataset
sns.lmplot(x='x', y='y', data=DF({'x':x,'y':y}))
plt.show()



[ alpha for r in res['raw_estimates'] for alpha in r['alphas'] if r['status'] != 'optimal' ]

sum( alpha['estimate'] for r in res['raw_estimates'] for alpha in r['alphas'] if r['alphas'] != 'optimal' )

res['raw_estimates']


res['basic_analysis'][0]["anion_approached_cation"]
res['intermediate_analysis'][0]["anion_approached_cation"]
res['advanced_analysis'][0]["anion_approached_cation"]



probs, counts = res['basic_analysis']


[ alpha for r in res['raw_estimates'] for alpha in r['alphas'] if alpha['q'] == 9 ]



counts['unreacted_precursors']
counts['unreacted_precursors']

denominator = sum( alpha['estimate'] for r in res['raw_estimates'] for alpha in r['alphas'] )
nominator = [ alpha['estimate'] for r in res['raw_estimates'] for alpha in r['alphas'] if alpha['q'] == 9 ][0]

1. - nominator/denominator




[ alpha for r in res['raw_estimates'] for alpha in r['alphas'] if alpha['molType'] == 'precursor']


mz, intensity = mol['spectrum']
DF({'mz':mz, 'intensity':intensity}).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/debugging/'+file_name +'csv',
    index = False
)

x = make_formulas('MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', Q = 1)

list(x.makeMolecules())[0][1]






short_data          = res['for_plot']['G_nodes_data']
remaining_peaks     = res['for_plot']['remaining_peaks']
long_data           = res['for_plot']['MIG_paths_data']

DF(short_data).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/test/short.csv',
    index = False )

DF(remaining_peaks).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/test/remaining_peaks.csv',
    index = False )

DF(long_data).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/test/long.csv',
    index = False )
