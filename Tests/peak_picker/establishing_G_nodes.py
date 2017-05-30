from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
from collections import Counter
import pandas as pd
import numpy as np


mol = substanceP.copy()
cut_off = 10.0; opt_P = 1.0

max_times_solve=30
joint_probability_of_envelope =.999
mz_prec=.05; prec_digits = 2
min_prob_of_envelope_in_picking = .7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursor_charge= mol['Q'],
                prec_digits     = prec_digits,
                joint_probability_of_envelope = joint_probability_of_envelope,
                mz_prec         = mz_prec,
                modifications   = mol['modifications']  )

# M.read_n_preprocess_spectrum(   spectrum    = mol['spectrum'],
#                                 cut_off     = 100.0 )

M.read_n_preprocess_spectrum(   spectrum    = mol['spectrum'],
                                opt_P       = 1.00     )

M.prepare_problems(min_prob_of_envelope_in_picking)


M.run(  solver  = 'sequential',
        method  = 'MSE',
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
        verbose = verbose )

SFG = M.res[0]['small_graph']
SFG.nodes(data=True)

SFG.edges(data=True)


vis_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/'

pd.DataFrame(M.spectra_iter()).to_csv(
        path_or_buf = vis_path+'spec_data.csv',
        index = False )

pd.DataFrame(M.i_get_original_spectrum()).to_csv(
        path_or_buf = vis_path+'original_spec.csv',
        index = False )


mz, intensity = M.spectrum
mz_trimmed, intensity_trimmed = M.spectrum_trimmed
mz = np.append(mz, mz_trimmed)
intensity = np.append(intensity, )

M.res[0]['small_graph'].nodes(data=True)

X = M.gen_ETDetective_inputs()
sum(X[1][k] for k in X[1])

print M.analyze_reactions('basic')
print M.analyze_reactions('inter')
print M.analyze_reactions('up_inter')


M.res



# from MassTodonPy.PeakPicker import PeakPicker
#
# PP = PeakPicker(M.Forms, M.IsoCalc)
# Graph_gen = PP.get_problems(M.spectrum)
# Graphs = list(Graph_gen)
#
# Graphs[0].nodes(data=True)
#
# from MassTodonPy.Deconvolutor import Deconvolutor_Min_Sum_Squares
#
#
# mol['spectrum'][1].sum()
# M.spectrum_intensity_trimmed
# M.spectrum_intensity
#
# for m in mol['spectrum'][1]:
#     print m
#
# DMSE = Deconvolutor_Min_Sum_Squares(Graphs[0])
# DMSE.run()
