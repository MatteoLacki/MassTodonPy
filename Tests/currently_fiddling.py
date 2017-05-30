from MassTodonPy import MassTodon, MassTodonize
from MassTodonPy.TestScripts import substanceP, ubiquitin
from collections import Counter
import pandas as pd
import numpy as np

#TODO: make simple substance P spectra comparison (fitting test)
#TODO: add info on the amounts of stuff making up the total G intensity
#TODO: run bootstrap
#TODO: run in sillico
#TODO: get threshold intensity


mol = substanceP.copy()
cut_off = 10.0
opt_P   = 1.0

max_times_solve = 30
joint_probability_of_envelope =.999
mz_prec =.05
min_prob_of_envelope_in_picking = .7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver  = 'sequential'
method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursor_charge= mol['Q'],
                joint_probability_of_envelope = joint_probability_of_envelope,
                mz_prec         = mz_prec,
                modifications   = mol['modifications']  )

M.read_n_preprocess_spectrum(
    spectrum    = mol['spectrum'],
    opt_P       = 1.00              )

M.spectra.keys()

M.prepare_problems(min_prob_of_envelope_in_picking)

M.run(  solver  = 'sequential',
        method  = 'MSE',
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
        verbose = verbose )


M.summarize_results()

M.res[0]['SG'].nodes(data=True)


M.analyze_reactions('basic')
M.analyze_reactions('inter')
M.analyze_reactions('up_inter')

pd.DataFrame(M.export_information_for_spectrum_plotting())
pd.DataFrame(M.export_information_for_spectrum_plotting(True))
