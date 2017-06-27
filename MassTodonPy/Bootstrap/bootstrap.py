from    numpy.random            import multinomial
from    itertools               import izip, product, repeat
from    MassTodonPy.PeakPicker  import PeakPickerBootstrap
from    multiprocessing         import Pool
from    MassTodonPy.Misc        import cvxopt_wrapper
import  numpy                   as np



def bootstrap_worker(args):
    '''This worker performs one whole bootstrap task.'''

    cut_off, original_total_intensity, mzs, intensities, mz_prec, clusters, ions_no, min_prob_per_molecule, Q, fasta, verbose = args

    if verbose:
        print 'Bootstrap worker online!'

    spectra = {}
    spectra['cut_off'] = cut_off
    spectra['original total intensity'] = original_total_intensity
    spectra_trimmed, spectra_left_out   = trim_spectrum(mzs, intensities, cut_off)
    spectra['total intensity after trim'] = spectra_trimmed[1].sum()
    spectra['trimmed intensity'] = original_total_intensity - spectra['total intensity after trim']
    spectrum_boot = dict(izip(*spectra_trimmed))
    PeakPickBoot  = PeakPickerBootstrap(mz_prec)
    problems = PeakPickBoot.turn_clusters_2_problems(clusters, spectrum_boot, ions_no, min_prob_per_molecule)

    res = solve(problems= problems,
                args    = { 'verbose' : verbose },
                solver  = 'sequential',
                multiprocesses_No = None,
                method  = 'MSE',
                max_times_solve = 5,
                verbose = verbose   )

    spectra['intensity of peaks paired with isotopologues'] = PeakPickBoot.stats['total intensity of experimental peaks paired with isotopologues']

    results = {}
    results['summary'] = summarize_results( spectra = spectra,
                                            raw_masstodon_res = res  )

    results.update(dict(map(lambda analyzer:
            (   analyzer,
                match_cz_ions(  results_to_pair         = res,
                                Q                       = Q,
                                fasta                   = fasta,
                                advanced_args           = {},
                                min_acceptEstimIntensity= min_prob_per_molecule,
                                analyzer                = analyzer,
                                accept_nonOptimalDeconv = True,
                                verbose                 = verbose  )
            ),['basic', 'intermediate', 'advanced'] )))

    if verbose:
        print 'Bootstrap worker finished!'

    return results


def run_bootstrap(  bootstrap_repeats,
                    spectra,
                    clusters,
                    mz_prec,
                    Q,
                    fasta,
                    min_prob_per_molecule = .75,
                    ions_no             = 100000,
                    multiprocesses_No   = None,
                    verbose             = False
    ):
    '''Run a bootstrap procedure for an existing MassTodon task.'''

    if verbose:
        print 'Running bootstrap!'
        print
    mzs, intensities = spectra['trimmed']
    boot_intensities = list(multinomial(ions_no, intensities/intensities.sum(), bootstrap_repeats).astype(np.float) * intensities.sum()/float(ions_no))

    cut_off = spectra['cut_off']
    original_total_intensity = spectra['original total intensity']

    iter_of_args = product(
        repeat(cut_off),
        repeat(original_total_intensity),
        repeat(mzs),
        boot_intensities,
        repeat(mz_prec),
        repeat(clusters), #TODO check if it works better to make out of this a global variable and TODO check if using Process directly works better.
        repeat(ions_no),
        repeat(min_prob_per_molecule),
        repeat(Q),
        repeat(fasta),
        repeat(verbose)
    )

    if verbose:
        print 'Problems ready'
        print

    with cvxopt_wrapper():
        P = Pool(multiprocesses_No)
        results = P.map( bootstrap_worker, iter_of_args )
        P.close()
        P.join()

    return results
