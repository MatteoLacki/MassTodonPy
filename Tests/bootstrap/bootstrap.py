import  pandas      as pd
import  cPickle     as pickle
import  json
from    bootstrap_misc import MassTodon_bootstrap
from    MassTodonPy     import MassTodonize, MassTodon
from    numpy.random    import multinomial
import  numpy           as np
from    itertools       import izip, repeat
from    multiprocessing import Pool

data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

ions_no         = 10**6
bootstrap_size  = 10

ID, mol = 0, substancesP[0]
mz_prec = .05
verbose = False
opt_P   = .99
# cut_off = 100
min_prob_of_envelope_in_picking = .7
method  = 'MSE'
solver  = 'sequential'
max_times_solve = 10
L1_x = L2_x = L1_alpha = L2_alpha = .001
multiprocesses_No = None
min_prob_per_molecule = .75

WH = mol['experimental_setting']['WH']
WV = mol['experimental_setting']['WV']
spectrum = mol['spectrum']

M = MassTodon(  fasta           = mol['fasta'],
                precursor_charge= mol['precursorCharge'],
                mz_prec         = mz_prec,
                modifications   = mol['modifications'],
                verbose         = verbose   )

M.read_n_preprocess_spectrum(   spectrum = spectrum,
                                opt_P    = opt_P    )

######### This will proceed with the orginal calculations.
M.run(  solver  = solver,
        multiprocesses_No = multiprocesses_No,
        method  = method,
        max_times_solve = max_times_solve,
        min_prob_per_molecule = min_prob_per_molecule,
        bootstrap = True,
        L1_x = L1_x,
        L2_x = L2_x,
        L1_alpha = L1_alpha,
        L2_alpha = L2_alpha,
        verbose  = verbose    )

Results = {}
Results['summary']              = M.summarize_results()
Results['basic analysis']       = M.analyze_reactions('basic')
Results['intermediate analysis']= M.analyze_reactions('intermediate')
Results['advanced analysis']    = M.analyze_reactions('advanced')
# Results['summary']

######### Here we will do bootstrap (actually we could add original tasks to the pool.. the hell with it...)
mzs, intensities = M.spectra['trimmed']

boot_intensities = multinomial(ions_no, intensities/intensities.sum(), bootstrap_size).astype(np.float) * intensities.sum()/float(ions_no)
clusters = M.clusters

from MassTodonPy.PeakPicker import PeakPickerBootstrap
from MassTodonPy.Solver     import solve
from MassTodonPy.Parsers    import trim_spectrum
from MassTodonPy.Summarator import summarize_results
from MassTodonPy.MatchMaker import match_cz_ions
from MassTodonPy.Misc       import cvxopt_wrapper



args = M.spectra['cut_off'], M.spectra['original total intensity'], (mzs, boot_intensities[0,]), mz_prec, clusters, ions_no, min_prob_per_molecule, mol['precursorCharge'], mol['fasta'], verbose

def bootstrap_worker(args):
    '''This worker performs one whole bootstrap task.'''
    cut_off, original_total_intensity, spectrum, mz_prec, clusters, ions_no, min_prob_per_molecule, Q, fasta, verbose = args
    mzs, intensities = spectrum

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

    Results = {}
    Results['summary'] = summarize_results( spectra = spectra,
                                            raw_masstodon_res = res  )

    Results.update(dict(map(lambda analyzer:
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
    return Results

P = Pool(multiprocesses_No)
results = P.map( bootstrap_worker, xrange(bootstrap_size) )
P.close()
P.join()

print len(results)
