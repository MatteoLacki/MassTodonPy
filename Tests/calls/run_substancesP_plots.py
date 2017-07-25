import cPickle as pickle
import os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy                            import MassTodonize
from MassTodonPy.TestScripts.substancesP    import substancesP
from MassTodonPy.Outputing.to_etdetective   import results_to_etdetective
from time                                   import time
from pandas                                 import DataFrame as DF


output_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/substancesP/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

# mol = substancesP[0]; ID = 0
for ID, mol in enumerate(substancesP):
    fasta= mol['fasta']
    Q    = mol['precursorCharge']
    res  = MassTodonize(fasta           = fasta,
                        precursor_charge= Q,
                        mz_prec         = .1,
                        joint_probability_of_envelope = .999,
                        spectrum        = mol['spectrum'],
                        modifications   = mol['modifications'],
                        opt_P           = .99,
                        solver          = 'multiprocessing',
                        multiprocesses_No = None,
                        max_times_solve = 10,
                        raw_data        = False,
                        for_plot        = True,
                        verbose         = True )

    WH, WV = [ mol['experimental_setting'][x] for x in ('WH','WV') ]
    file_path = output_path + str(ID) + '_WH-' + str(WH) + '_WV-' + str(WV) + '/'

    short_data          = res['for_plot']['G_nodes_data']
    remaining_peaks     = res['for_plot']['remaining_peaks']
    long_data           = res['for_plot']['MIG_paths_data']

    if not os.path.exists(file_path):
        os.makedirs(file_path)

    DF(short_data).to_csv(  path_or_buf = file_path + 'short.csv',
                            index = False )

    DF(remaining_peaks).to_csv( path_or_buf = file_path + 'remaining_peaks.csv',
                                index = False )

    DF(long_data).to_csv(   path_or_buf = file_path + 'long.csv',
                            index = False )
