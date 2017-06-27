import  cPickle     as pickle
from    itertools       import izip, repeat, islice
from    multiprocessing import Pool
from    numpy.random    import multinomial
import  numpy as np
from    MassTodonPy import MassTodonize


def bootstrap_worker(worker_args):
    info, fasta, Q, mz_prec, opt_P, modifications, masses, intensities, max_times_solve, verbose = worker_args
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
    results.update(info)
    return results


def bootstrap_substance_P(  mol,
                            bootstrap_size  = 1000,
                            ions_no         = 10**6,
                            mz_prec         = 0.05,
                            opt_P           = .999,
                            max_times_solve = 10,
                            multiprocesses_No = None,
                            verbose         = False
    ):
    '''Run a bootstrap.'''
    mzs, intensities = mol['spectrum']
    sim_intensities  = multinomial(ions_no, intensities/intensities.sum(), bootstrap_size).astype(np.float) * intensities.sum()/float(ions_no)
    pool_args = izip(   repeat(mol['info']),
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


def analyze_experiments(substances,
                        results_path,
                        K = None,
                        bootstrap_size  = 1000,
                        ions_no         = 10**6,
                        mz_prec         = 0.05,
                        opt_P           = .999,
                        max_times_solve = 10,
                        multiprocesses_No = None,
                        verbose         = False  ):
    for ID, mol in enumerate(islice(substances, K)):
        results = {}
        WH = mol['experimental_setting']['WH']
        WV = mol['experimental_setting']['WV']
        mol['info'] = {'WH':WH, 'WV':WV, 'ID':ID}

        results['real'] = bootstrap_worker((mol['info'],
                                            mol['fasta'],
                                            mol['precursorCharge'],
                                            mz_prec,
                                            opt_P,
                                            mol['modifications'],
                                            mol['spectrum'][0],
                                            mol['spectrum'][1],
                                            max_times_solve,
                                            verbose))

        results['bootstrap'] = bootstrap_substance_P(mol, bootstrap_size)
        with open(results_path+str(ID), 'w') as handler:
            pickle.dump(results, handler)
        print 'Dumped', ID, 'out of', len(substances)
