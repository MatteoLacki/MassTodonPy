from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
from pandas import DataFrame as DF
import cPickle as pickle

# with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubiquitins.example', 'r') as h:
#     ubiquitins = pickle.load(h)
# mol = ubiquitins[120]

mol = substanceP.copy()

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    raw_data        = True,
                    for_plot        = True,
                    verbose         = False )

short_data          = res['for_plot']['G_nodes_data']
remaining_peaks     = res['for_plot']['remaining_peaks']
long_data           = res['for_plot']['MIG_paths_data']

DF(short_data).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_short.csv',
    index = False )

DF(remaining_peaks).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tsests/visual/ubi_remaining_peaks.csv',
    index = False )

DF(long_data).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_long.csv',
    index = False )

# res['spectra'].keys()
# res['spectra']['']
# SG = res['raw_estimates'][0]['SG']
# SG.nodes(data=True)
#
# SG['I724']
# SG.nodes(data=True)[10]
# [g for g in SG.nodes(data=True) if g[0][0] == 'M']
# SG.edges(data=True)
