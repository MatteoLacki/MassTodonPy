%load_ext autoreload
%load_ext line_profiler
%autoreload
import  networkx as nx
from    Visualization import plot_spectrum, plot_deconvolution_graph
import  matplotlib.pyplot as plt
from    cvxopt import matrix, spmatrix, sparse, spdiag, solvers
from    PeakPicker import create_G_nodes, getGraphs, trim_unlikely_molecules, contains_experimental_peaks, add_zero_intensity_G_nodes
import  cPickle as pickle
from    Solver import prepare_deconvolution

G = nx.Graph()
G.add_node('M1', type='M', formula='C52H89N13O13S1', g=-1, molType='c10', q=1)
G.add_node('M2', type='M', formula='C52H89N13O13S1', g=2,  molType='c12', q=3)
G.add_node('I11', type='I', intensity=.4)
G.add_node('I12', type='I', intensity=.4)
G.add_node('I13', type='I', intensity=.2)
G.add_node('I21', type='I', intensity=.6)
G.add_node('I22', type='I', intensity=.2)
G.add_node('I23', type='I', intensity=.2)
G.add_node('G1', type='G', intensity=1000.)
G.add_node('G2', type='G', intensity=600.)
G.add_node('G3', type='G', intensity=200.)
G.add_node('G4', type='G', intensity=200.)
G.add_edges_from(
    [   ('M1','I11'),('M1','I12'),('M1','I13'),
        ('M2','I21'),('M2','I22'),('M2','I23'),
        ('I11','G1'),('I12','G2'),('I13','G3'),
        ('I21','G1'), ('I22','G2'),('I23','G4')])
add_zero_intensity_G_nodes(G)
# plot_deconvolution_graph(G)

L2_percent = 0.000001

SFG = G.copy()
total_G_intensity, squared_G_intensity, cnts, varNo, P_mat, q_vec, G_mat, h_vec, A_mat, b_vec, initvals = prepare_deconvolution(SFG, L2_percent)
cnts
varNo
total_G_intensity
squared_G_intensity
# print P_mat
# , q_vec, G_mat, h_vec, A_mat, b_vec, initvals


solvers.options['show_progress'] = False
sol = solvers.qp(P_mat,q_vec,G_mat,h_vec,A_mat,b_vec, initvals=initvals)


alphas = []
for N in SFG:
    type_N = SFG.node[N]['type']
    if type_N == 'M':
        estimated_alpha = sol['x'][ cnts['GI']+SFG.node[N]['cnt'] ]
        SFG.node[N]['estimate'] = estimated_alpha
        alphas.append( SFG.node[N].copy() )
    if type_N == 'G':
        for I in SFG[N]:
            estimated_x = sol['x'][ SFG.edge[N][I]['cnt'] ]
            SFG.edge[N][I]['estimate'] = estimated_x

SFG.nodes(data=1)
SFG.edge
alphas

error = 0.0
for G in SFG:
    if SFG.node[G]['type'] == 'G':
        I_intensity = SFG.node[G]['intensity']
        outflow = 0.0
        for I in SFG[G]:
            outflow += sol['x'][ SFG.edge[G][I]['cnt'] ]
        error += (I_intensity - outflow)**2

from math import sqrt
error = sqrt(error)

error
