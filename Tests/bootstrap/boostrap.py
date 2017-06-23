import  pandas      as pd
import  cPickle     as pickle
import  json
from    bootstrap_misc import MassTodon_bootstrap
from    MassTodonPy     import MassTodonize, MassTodon
from    numpy.random    import multinomial
import  numpy as np
from    itertools import izip
from    multiprocessing import Pool

data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

ions_no         = 10**6
bootstrap_size  = 1

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
from MassTodonPy.PeakPicker import PeakPickerBootstrap

clusters = M.provide_clusters()
spectrum = M.spectra['trimmed']

PeakPickBoot = PeakPickerBootstrap(mz_prec)
problems = PeakPickBoot.turn_clusters_2_problems(clusters, spectrum, ions_no, min_prob_per_molecule)



from MassTodonPy.Solver import solve


solve(  problems = problems,
        args   = args,
        solver = solver,
        multiprocesses_No = multiprocesses_No,
        method = method,
        max_times_solve = max_times_solve,
        verbose= self.verbose   )






# multiprocesses_No = None
P = Pool(multiprocesses_No)
results = P.map(    helper,
                    pool_args    )

P.close()
P.join()
