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
from    PeakPicker import create_G_nodes, getGraphs, trim_unlikely_molecules, contains_experimental_peaks
import  cPickle as pickle
from    Solver import prepare_deconvolution
# import  multiprocessing as multiKulti

path    = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
spectrum= ParseMzXML(path+'Ubiquitin_ETD_10 ms_1071.mzXML',cut_off_intensity=100)
    # plot_spectrum(spectrum, 1215, 1230) # plot_spectrum(spectrum, 0,8000)
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8; modifications = {}; jP = .999; mzPrec = .05; precDigits = 2

massTodon = MassTodon(  fasta           = fasta,
                        precursorCharge = Q,
                        precDigits      = precDigits,
                        mzPrec          = mzPrec )

recalculate   = True
problems_file = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/problems.p"


if recalculate:
    BFG     = massTodon.peakPicker.BFG_representation(spectrum)
    css     = nx.connected_component_subgraphs(BFG)
    problems= list(getGraphs(css, 0.9))
    pickle.dump( problems, open( problems_file, "wb" ) )
problems = pickle.load( open( problems_file, "rb" ) )

for P in problems:
    create_G_nodes(P)

def deconvolve(SFG, L2_percent=0.0):
    SFG = SFG.copy()
    cnts, varNo, G_intensity, P_mat, q_vec, G_mat, h_vec, A_mat, b_vec, initvals = prepare_deconvolution(SFG, L2_percent)
    sol = solvers.qp(P_mat,q_vec,G_mat,h_vec,A_mat,b_vec, initvals=initvals)
    return sol

solvers.options['show_progress'] = False
%%time
sols = [deconvolve(P, 0.0) for P in  problems]

unknowns = [ P for s,P in zip(sols,problems) if s['status']=='unknown' ]
len(unknowns)
U = unknowns[-3]
# [len(u) for u in unknowns]
len(unknown)
plot_deconvolution_graph(unknown)


for G in U:
    if U.node[G]['type']=='G':
        print U.node[G]





def get_total_intensity(G):
    G = G.copy()
    add_missing_experimental_groups(G)
    return prepare_deconvolution(G)[0]

get_total_intensity(unknown)


totalI = [ get_total_intensity(P) for P in problems ]

solutions_n     = Counter()
total_intensity = Counter()
for intensity, sol in zip(totalI, sols):
    solutions_n[sol['status']] += 1
    total_intensity[sol['status']] += intensity

solutions_n

for t in solutions_n.keys():
    total_intensity[t] /= solutions_n[t]

total_intensity['optimal']/total_intensity['unknown']*100
Counter( s['status'] for s in sols)





# P = multiKulti.Pool(3)
# %%time
# sols = P.map( deconvolve, problems )

Counter( s['status'] for s in sols)
Counter( s['iterations'] for s in sols if s['status'] == 'unknown')
Counter( s['iterations'] for s in sols if s['status'] == 'optimal')

[s['x'] for s in sols if s['status'] == 'unknown']

[s['x'] for s in sols if s['status'] == 'optimal']


def check_info(SFG):
    SFG = SFG.copy()
    add_missing_experimental_groups(SFG)
    return prepare_deconvolution(SFG)

infos = [check_info(P) for P in problems]
info = infos[0]

from scipy.linalg import norm
P_mat = np.matrix(matrix(info[3]))
norm(P_mat,ord=2)


def check_rank(SFG):
    SFG = SFG.copy()
    add_missing_experimental_groups(SFG)
    cnts, varNo, G_intensity, P_mat, q_vec, G_mat, h_vec, A_mat, b_vec, initvals = prepare_deconvolution(SFG)
    rank = np.linalg.matrix_rank(np.array(matrix(A_mat)))
    size = A_mat.size
    return rank, size[0], size[1]

ranks = [check_rank(P) for s, P in zip(sols,problems) if s['status'] == 'unknown']
ranks







# ordering some of the problem graph nodes and edges
smallProblems = [P for P in problems if len(P) < 20]
doubleMolProblem = [P for P in problems if len([N for N in P if P.node[N]['type']=='M']) == 2]
tripleMolProblem = [P for P in problems if len([N for N in P if P.node[N]['type']=='M']) == 3]
# problem = problems[0].copy()
# problem = smallProblems[0].copy()
# problem = doubleMolProblem[0].copy()
problem = tripleMolProblem[0].copy()
add_missing_experimental_groups(problem)
# plot_deconvolution_graph(problem)
problem


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
        twos = matrix( 1.0, (G_degree,1))
        P_list.append( 2.0 * twos * twos.T ) # 2 because of .5 x'Px parametrization.

q_list.append(  matrix( 0.0, ( cnts['M'], 1 ) )  )
q_vec = matrix(q_list)
P_list.append(  matrix( 0.0, ( cnts['M'], cnts['M'] ) )  )
P_spectral_norm = 2.0 * max(p_mat.size[0] for  p_mat in P_list) # spec(11')=dim 1
P_mat = spdiag(P_list)/P_spectral_norm
q_vec = q_vec/P_spectral_norm
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

A_mat   = spmatrix( A_x, A_i, A_j, size=( cnts['I'], varNo ) )
b_vec   = matrix( 0.0, ( cnts['I'], 1)  )
initvals= {}
initvals['x'] = matrix( 0.0, ( varNo, 1)  )


def normalize_rows(M):
    for i in xrange(M.size[0]):
        row_hopefully = M[i,:]
        M[i,:] = row_hopefully/sum(row_hopefully)

normalize_rows(A_mat)




solvers.options['show_progress'] = True
sol = solvers.qp(P_mat,q_vec,G_mat,h_vec,A_mat,b_vec, initvals=initvals)
sol = solvers.qp(P_mat,q_vec,G_mat,h_vec,A_mat,b_vec)

print sol['x']

sol['status']

varNo

# solvers.options['show_progress'] = False

np.array(matrix(A_mat))


#
# cnts
B = np.array([[1,1,-.05,0,0,0],[0,0,0,1,1,-.4]])
# np.linalg.matrix_rank(B)
#
# BB = np.array([[1,1,-.05],[0,0,-.4]])
# print BB
# np.linalg.matrix_rank(BB)
