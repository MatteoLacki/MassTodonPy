%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon import MassTodon
import  numpy as np
from    math import sqrt
from    pandas import DataFrame
from    frozendict import frozendict
from    collections import Counter, defaultdict
import  networkx as nx
import  igraph as ig
from    Parsers import ParseMzXML
from    Visualization import plot_spectrum, plot_deconvolution_graph
import  matplotlib.pyplot as plt
from    cvxopt import matrix, spmatrix, sparse, spdiag, solvers
from    PeakPicker import group_experimental_peaks, getGraphs, trim_unlikely_molecules, contains_experimental_peaks
import cPickle as pickle

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
spectrum = ParseMzXML(path+'Ubiquitin_ETD_10 ms_1071.mzXML',cut_off_intensity=100)
# plot_spectrum(spectrum, 1215, 1230) # plot_spectrum(spectrum, 0,8000)
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8; modifications = {}; jointProb = .999; mzPrec = .05;
massTodon = MassTodon(fasta, Q, massPrecDigits=2)

# BFG     = massTodon.peakPicker.BFG_representation(spectrum)
# css     = nx.connected_component_subgraphs(BFG)
# problems= list(getGraphs(css, 0.7))

problems_file = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/problems.p"
# pickle.dump( problems, open( problems_file, "wb" ) )
problems = pickle.load( open( problems_file, "rb" ) )

for P in problems:
    group_experimental_peaks(P)

smallProblems = [P for P in problems if len(P) < 20]
# problem = problems[0].copy()
problem = smallProblems[0].copy()

    # adding missing G edges with 0 intensity
cnts = Counter(problem.node[N]['type'] for N in problem)

Gcnt = cnts['G']
newGnodes = []
newGIedges= []
for I in problem:
    if problem.node[I]['type'] == 'I':
        if len( problem[I] ) == 1:
            G = 'G' + str(Gcnt)
            newGnodes.append( (G,{'intensity':0.0, 'type':'G'}) )
            newGIedges.append( (G,I) )
            Gcnt += 1

problem.add_nodes_from(newGnodes)
problem.add_edges_from(newGIedges)

# ordering some of the problem graph nodes and edges
cnts = Counter()
for N in problem:
    Ntype = problem.node[N]['type']
    problem.node[N]['cnt'] = cnts[Ntype]
    cnts[Ntype] += 1
    if Ntype == 'G':
        for I in problem.edge[N]:
            problem.edge[N][I]['cnt'] = cnts['GI']
            cnts['GI'] += 1

varNo = cnts['GI']+cnts['M']

squared_G_intensity = 0.0
q_list = []
P_list = []
for G in problem:
    if problem.node[G]['type']=='G':
        G_intensity = problem.node[G]['intensity']
        squared_G_intensity += G_intensity
        G_degree = len(problem[G])
        q_list.append(  matrix( -G_intensity, size=(G_degree,1) )  )
        twos = matrix(2.0, (G_degree,1)) # 2 because of .5 x'Px parametrization.
        P_list.append( twos * twos.T )

q_list.append(  matrix( 0.0, ( cnts['M'], 1 ) )  )
q_vec = matrix(q_list)
P_list.append(  matrix( 0.0, ( cnts['M'], cnts['M'] ) )  )
P_mat = spdiag(P_list)
G_mat = spmatrix(-1.0, xrange(varNo), xrange(varNo))
h_vec = matrix(0.0, size=(varNo,1) )

A_x = [];A_i = [];A_j = []
probabilities = Counter()
isotopologueNo= Counter()
for M in problem:
    if problem.node[M]['type']=='M':
        M_cnt = problem.node[M]['cnt']
        for I in problem[M]:
            i_cnt = problem.node[I]['cnt']
            probabilities[M] += problem.node[I]['intensity']
            A_x.append(-problem.node[I]['intensity'] )
            A_i.append( i_cnt )
            A_j.append( M_cnt + cnts['GI'] )
            isotopologueNo[M_cnt + cnts['GI']] += 1
            for G in problem[I]:
                if not G == M:
                    A_x.append( 1.0 )
                    A_i.append( i_cnt )
                    A_j.append( problem.edge[G][I]['cnt'] )

A_mat = spmatrix( A_x, A_i, A_j, size=( cnts['I'], varNo ) )
b_vec = matrix( 0.0, ( cnts['I'], 1)  )
sol = solvers.qp(P_mat,q_vec,G_mat,h_vec,A_mat,b_vec)


A_mat.size
print A_mat

plot_deconvolution_graph(problem)

# A
# A_np = np.matrix(matrix(A))
# np.linalg.matrix_rank(A_np)
#
# cnts
# B = np.array([[1,1,-.05,0,0,0],[0,0,0,1,1,-.4]])
# np.linalg.matrix_rank(B)
#
# BB = np.array([[1,1,-.05],[0,0,-.4]])
# print BB
# np.linalg.matrix_rank(BB)
