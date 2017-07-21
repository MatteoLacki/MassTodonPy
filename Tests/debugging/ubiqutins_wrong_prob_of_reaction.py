from MassTodonPy import MassTodonize
import cPickle as pickle
with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)

# mol = ubiquitins[120]

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    highcharts      = True,
                    verbose         = False )
