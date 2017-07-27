import  os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy import MassTodonize
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
import cPickle as pickle
from pandas import DataFrame as DF

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)

output_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiquitins/'

# ubiquitin = ubiquitins[0]
for ubiquitin in ubiquitins:
    fasta= ubiquitin['fasta']
    Q    = ubiquitin['Q']
    res  = MassTodonize(fasta           = fasta,
                        precursor_charge= Q,
                        mz_prec         = .05,
                        joint_probability_of_envelope = .999,
                        spectrum        = ubiquitin['spectrum'],
                        opt_P           = .95,
                        solver          = 'multiprocessing',
                        multiprocesses_No = None,
                        max_times_solve = 10,
                        raw_data        = False,
                        for_plot        = True,
                        verbose         = True )


    file_path = output_path + ubiquitin['experimental_setting']['files'].split('.')[0] + '/'

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


# from MassTodonPy.Formulator import make_formulas
# list(make_formulas(ubiquitin['fasta'], ubiquitin['Q']).makeMolecules())
