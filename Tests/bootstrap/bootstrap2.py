import  pandas      as pd
import  cPickle     as pickle
import  json
from    itertools       import izip, repeat, islice
from    multiprocessing import Pool
from    numpy.random    import multinomial
import  numpy as np
from    MassTodonPy import MassTodonize


def bootstrap_worker(worker_args):
    WH, WV, ID, fasta, Q, mz_prec, opt_P, modifications, masses, intensities, max_times_solve, verbose = worker_args
    results = MassTodonize(
        fasta           = fasta,
        precursor_charge= Q,
        mz_prec         = mz_prec,
        opt_P           = opt_P,
        modifications   = modifications,
        spectrum        = (masses, intensities),
        solver          = 'sequential',
        max_times_solve = max_times_solve,
        raw_data        = False,
        verbose         = verbose )
    results['WH'] = WH
    results['WV'] = WV
    results['ID'] = ID
    return results


def bootstrap(  mol,
                ID,
                bootstrap_size  = 1000,
                ions_no         = 10**6,
                mz_prec         = 0.05,
                opt_P           = .999,
                max_times_solve = 10,
                multiprocesses_No = None,
                verbose         = False
    ):
    '''Run a bootstrap.'''
    WH = mol['experimental_setting']['WH']
    WV = mol['experimental_setting']['WV']
    mzs, intensities = mol['spectrum']
    sim_intensities = multinomial(ions_no, intensities/intensities.sum(), bootstrap_size).astype(np.float) * intensities.sum()/float(ions_no)
    pool_args = izip(   repeat(WH),
                        repeat(WV),
                        repeat(ID),
                        repeat(mol['fasta']),
                        repeat(mol['precursorCharge']),
                        repeat(mz_prec),
                        repeat(opt_P),
                        repeat(mol['modifications']),
                        repeat(mzs),
                        sim_intensities,
                        repeat(max_times_solve),
                        repeat(verbose)             )
    P = Pool(multiprocesses_No)
    results = P.map( bootstrap_worker, pool_args )
    P.close()
    P.join()
    return results


data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'

results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_CSV_27_06_2017/'

with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

K = 2
bootstrap_size = 10
# K = len(substancesP)

for ID, mol in enumerate(islice(substancesP, K)):
    results = bootstrap(mol, ID, bootstrap_size)
    with open(results_path+str(ID), 'w') as handler:
        pickle.dump(results, handler)
    print 'Dumped', ID, 'out of',len(substancesP),'.'
