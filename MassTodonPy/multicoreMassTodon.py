# %load_ext autoreload
# %load_ext line_profiler
%autoreload

from    MassTodon import MassTodon
import  numpy as np
from    intervaltree import Interval as I, IntervalTree as Itree
from    math import sqrt
from    pandas import DataFrame
from    frozendict import frozendict
from    collections import Counter
import  networkx as nx
import  igraph as ig
from    Parsers import ParseMzXML
from    Visualization import plot_spectrum, plot_connected_component

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
spectrum = ParseMzXML(path+'Ubiquitin_ETD_10 ms_1071.mzXML', cut_off_intensity=10)

# plot_spectrum(spectrum, 1215, 1230)
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8
modifications = {}
jointProb = .999
mzPrec = .05

massTodon = MassTodon(fasta, Q, massPrecDigits=2)
massTodon.peakPicker.setMassSpectrum(spectrum)

# approach focussed on Divide ed Impera
%%time

ePeaks = Itree(I(mz - mzPrec, mz + mzPrec, (mz, intensity)) for
               mz, intensity in np.rollaxis(massTodon.peakPicker.MS, 0))

edges = []
nodes = []
iso_cnt = 0
for cnt, (mType, formula, aaNo, q, g) in enumerate(massTodon.formulator.makeMolecules()):
    iso_mzs, iso_intensities = massTodon.isoCalc.isoEnvelope(formula, jointProb, q, g)
    M = mType + '_' + str(q) + '_' + str(g)
    nodes.append(())
    for isoMZ, isoI in zip(iso_mzs, iso_intensities):
        I = 'I'+str(iso_cnt)
        # nodes.append({ 'name':I, 'mz':isoMZ, 'intensity':isoI })
        edges.append((M,I))
        for exp in ePeaks[isoMZ]:
            expMZ, expI = exp.data
            edges.append((I,expMZ))
            exp_cnt += 1
        iso_cnt += 1

BFG = nx.Graph(edges)
ccs = list(nx.connected_component_subgraphs(G))

len(ccs)


# %%time
# edges = []
# nodes = []
# iso_cnt = 0
# exp_cnt = 0
# for cnt, (mType, formula, aaNo, q, g) in enumerate(massTodon.formulator.makeMolecules()):
#     iso_mzs, iso_intensities = massTodon.isoCalc.isoEnvelope(formula, jointProb, q, g)
#     missed_intensity = 0.0
#     M = mType + '_' + str(q) + '_' + str(g)
#     nodes.append({'name':M})
#     for isoMZ, isoI in zip(iso_mzs, iso_intensities):
#         I = 'I'+str(iso_cnt)
#         nodes.append({ 'name':I, 'mz':isoMZ, 'intensity':isoI })
#         edges.append({ 'source':M, 'target':I })
#         for exp in ePeaks[isoMZ]:
#             expMZ, expI = exp.data
#             E = 'E'+str(exp_cnt)
#             nodes.append({ 'name':E, 'mz':expMZ, 'intensity':expI })
#             exp_cnt += 1
#         iso_cnt += 1
#
# G = ig.Graph.DictList(vertices=nodes, edges=edges, iterative=False)
#
# clust = G.clusters()
# largest = G.clusters().giant()
# largest
# X = G.vs(clust)
# x = G.components().subgraphs()
# connected_components = [ cc for cc in G.decompose()]
