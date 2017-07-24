import  cPickle     as pickle
from    itertools       import izip, repeat, islice
from    multiprocessing import Pool
from    numpy.random    import multinomial
import  numpy as np
from    MassTodonPy import MassTodonize
import  traceback

def bootstrap_worker(worker_args):
    info, fasta, Q, mz_prec, opt_P, cut_off, modifications, masses, intensities, max_times_solve, verbose = worker_args
    try:
        results = MassTodonize(
            fasta           = fasta,
            precursor_charge= Q,
            mz_prec         = mz_prec,
            opt_P           = opt_P,
            cut_off         = cut_off,
            modifications   = modifications,
            spectrum        = (masses, intensities),
            solver          = 'sequential',
            max_times_solve = max_times_solve,
            raw_data        = False,
            verbose         = verbose )
        results.update(info)
    except Exception as e:
        traceback.print_exc()
        print()
        results = None
    return results


def bootstrap(  mol,
                bootstrap_size  = 1000,
                ions_no         = 10**6,
                mz_prec         = 0.05,
                opt_P           = None,
                cut_off         = None,
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
                        repeat(cut_off),
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
                        K               = None,
                        bootstrap_size  = 1000,
                        ions_no         = 10**6,
                        mz_prec         = 0.05,
                        cut_off         = None,
                        opt_P           = None,
                        max_times_solve = 10,
                        multiprocesses_No = None,
                        verbose         = False  ):

    for ID, mol in enumerate(islice(substances, K)):
        results  = {}
        settings = mol['experimental_setting']
        del mol['experimental_setting']
        settings['ID'] = ID
        mol['info'] = settings

        if verbose:
            print 'starting calculations on real data'
            print 'fasta:', mol['fasta']
            print 'Q:', mol['precursorCharge']

        results['real'] = bootstrap_worker((mol['info'],
                                            mol['fasta'],
                                            mol['precursorCharge'],
                                            mz_prec,
                                            opt_P,
                                            cut_off,
                                            mol['modifications'],
                                            mol['spectrum'][0],
                                            mol['spectrum'][1],
                                            max_times_solve,
                                            verbose))

        if verbose:
            print 'starting bootstrap calculations'

        results['bootstrap'] = bootstrap(
            mol,
            bootstrap_size,
            ions_no,
            mz_prec,
            opt_P,
            cut_off,
            max_times_solve,
            multiprocesses_No,
            verbose )

        with open(results_path+str(ID), 'w') as handler:
            pickle.dump(results, handler)
        print 'Dumped', ID+1, 'out of', len(substances)
