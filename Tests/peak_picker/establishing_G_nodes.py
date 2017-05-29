from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
from collections import Counter
import pandas as pd
import numpy as np


mol = substanceP.copy()
cut_off = 10.0; opt_P = 1.0

max_times_solve=30
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursorCharge = mol['Q'],
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = mol['modifications']  )

M.read_spectrum(spectrum = mol['spectrum'])

# mz, intensity = mol['spectrum']
#
#
#
# def trim_spectrum(mz, intensity, cutOff=cutOff):
#     '''Remove peaks below a given cut off.'''
#     mz_trimmed = mz[intensity <= cutOff]
#     mz = mz[intensity > cutOff]
#     intensity_trimmed = intensity[intensity <= cutOff]
#     intensity = intensity[intensity > cutOff]
#     return (mz, intensity) , (mz_trimmed, intensity_trimmed)
#
# trim_spectrum(mz, intensity, cutOff=100)

#
# Exps = Itree( II( mz-mzPrec, mz+mzPrec, (mz, intensity) ) for mz, intensity in zip(*massSpectrum) )

# list(M.Forms.makeMolecules())
# mol['spectrum'][1].max()
# print M.total_I
# print M.total_I_after_cut_off
# print M.total_I_after_percent_trim


M.prepare_problems(M_minProb)

# M.run(  solver  = 'sequential',
#         method  = 'MSE',
#         max_times_solve = max_times_solve,
#         L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
#         verbose = verbose )

M.run(  solver  = 'sequential',
        method  = 'MSE',
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
        verbose = verbose )

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
