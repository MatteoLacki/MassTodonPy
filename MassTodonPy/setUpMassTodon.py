%load_ext autoreload
%autoreload
%load_ext line_profiler

from MassTodon import MassTodon
import numpy as np
from Solver import solveProblem, isDeconvoProb

from Parsers import ParseMzXML

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
spectrum_file = path+'Ubiquitin_ETD_10 ms_1071.mzXML'
spectrum = ParseMzXML(spectrum_file)
spectrum

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; P = .999

massTodon = MassTodon(fasta, Q, massPrecDigits = 1)

masses, intensities, noise_masses, noise_intensities = massTodon.randomSpectrum(ionsNo)
mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)
massSpectrum = list(zip(mz,ints))

massTodon.peakPicker.setMassSpectrum(massSpectrum)
massTodon.peakPicker.add_M_n_E()

#TODO what goes wrong with graph generation: why so slow?
massTodon.peakPicker.add_I()
massTodon.peakPicker.add_eG()
G = massTodon.peakPicker.G

import networkx as nx

%%time
deconvs = [isDeconvoProb(g) for g in nx.connected_component_subgraphs(G)]

len(deconvs)
from collections import Counter
Counter(deconvs)

%%time
solutions = []
wrongGraphs = []
for g in nx.connected_component_subgraphs(G):
    if isDeconvoProb(g):
        try:
            solutions.append(solveProblem(g))
        except:
            wrongGraphs.append(g)

wrongGraphs
solutions[0]
len(solutions)



Counter(sol[0].success for sol in solutions)

poorSolutions = [sol for sol in solutions if not sol[0].success]

len(poorSolutions[0][3]['A_eq'])

np.array(poorSolutions[0][3]['A_eq']).shape
from itertools import chain

[ np.array(sol[3]['A_eq']).shape for sol in solutions ]

[ np.array(sol[3]['A_eq']).shape for sol in poorSolutions ]

#TODO
