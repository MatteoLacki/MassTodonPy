%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import  MassTodon
from    Formulator      import  makeFormulas
import  cPickle         as      pickle
import  networkx        as      nx

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# modifications = {   ('N',2) :       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }
# Forms = makeFormulas(fasta=fasta, Q=Q, fragType='cz')
# M = MassTodon(  fasta           = fasta,
#                 precursorCharge = Q,
#                 precDigits      = precDigits,
#                 jointProbability= jP,
#                 mzPrec          = mzPrec )
# path  = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
# path += 'Ubiquitin_ETD_10 ms_1071.mzXML'
# cutOff = 100; topPercent = .999
# M.readSpectrum(path=path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
# M.prepare_problems(M_minProb)
# mu=1e-5; lam=0.0; nu=0.001
# %%time
# res = M.run(solver='sequential', method='MSE', mu=mu, lam=lam, nu=0.001)
# with open(file_path, 'w') as f:
    # pickle.dump(res, f)
with open(file_path, 'rb') as f:
    res = pickle.load(f)

from collections import defaultdict, Counter


# Counter(x['molType'] for e,E,s in res for x in e)
res

def update_nators(mol, nominator=0.0, denominator=0.0):
    ETnoDs  = mol['g']
    PTRs    = Q-mol['q']-mol['g']
    nominator   += PTRs*mol['estimate']
    denominator += (ETnoDs+PTRs)*mol['estimate']
    return nominator, denominator

import numpy as np
import scipy as spy
from scipy.optimize import linprog as SimplexAlgorithm

# denominator_PTR = nominator_PTR = 0.0
# fragments = defaultdict(lambda:defaultdict(Counter))
# L = len(fasta)
#
# for mols, error, status in res:
#     if status=='optimal':
#         for mol in mols:
#             if mol['molType']=='precursor':
#                 if mol['q']==Q and mol['g']==0:
#                     precursor_no_reaction = mol['estimate'] # this might be missing. use a dictionary and a absent key exception
#                 else:
#                     nominator_PTR, denominator_PTR = update_nators(mol, nominator, denominator)
#             else:
#                 mt = mol['molType']
#                 if mt[0]=='c':
#                     breakpoint=int(mt[1:])
#                 if mt[0]=='z':
#                     breakpoint=L-int(mt[1:])
#                 q=mol['q']
#                 fragments[breakpoint][mt[0]][(mt, q)]+= mol['estimate']

denominator_PTR = nominator_PTR = 0.0
FG = nx.Graph()
L = len(fasta)

for mols, error, status in res:
    if status=='optimal':
        for mol in mols:
            if mol['molType']=='precursor':
                if mol['q']==Q and mol['g']==0:
                    precursor_no_reaction = mol['estimate'] # this might be missing. use a dictionary and a absent key exception
                else:
                    nominator_PTR, denominator_PTR = update_nators(mol, nominator, denominator)
            else:
                frag = (mol['molType'],mol['q'])
                FG.add_node(frag, intensity=mol['estimate'] )
                FG.add_edge(frag,frag)

prob_PTR = nominator_PTR/denominator_PTR
prob_ETnoD = 1.0 - prob_PTR

for C, qC in FG:
    if C[0]=='c':
        for Z, qZ in FG:
            if Z[0]=='z':
                bpC = int(C[1:])
                bpZ = L - int(Z[1:])
                if bpC==bpZ and qC + qZ < Q-1:
                    FG.add_edge((C,qC),(Z,qZ))

import  pylab as pl
from    matplotlib import collections  as mc
import  matplotlib.pyplot as plt
FG

nx.draw(FG, node_size=10, pos=nx.random_layout(FG))
plt.show()

ccs = [cc for cc in nx.connected_component_subgraphs(FG)]

[ (i,len(cc)) for i,cc in enumerate(ccs) ]

ccs[23].nodes()
ccs[22].nodes()
ccs[23].edges()

frags = fragments[2]

ccs[22].nodes(data=True)
frags.node[('z45',3)]['intensity']

frags.nodes()
frags = ccs[0]

[collect_fragments(cc,Q) for cc in ccs]


frags = ccs[23]

frags.edges()


# numbering edges
for i, E in enumerate(frags.edges(data=True)):
    frags.edge[E[0]][E[1]]['cnt'] = i

# finding costs
c = []
for A, B in frags.edges():
    if A == B:
        c.append(Q-1-A[1])
    else:
        c.append(Q-1-A[1]-B[1])

c = np.array(c)
b = np.array([ x[1]['intensity'] for x in frags.nodes(data=True) ])

def get_incidence_matrix(frags):
    A = np.zeros((len(frags), len(frags.edges())))
    for i, N in enumerate(frags):
        for M in frags.edge[N]:
            j = frags.edge[N][M]['cnt']
            A[i,j] = 1.0
    return A

A = get_incidence_matrix(frags)

A, b, c

optim_result = SimplexAlgorithm(c, A_eq=A, b_eq=b, options={"disp": True})




# def collect_fragments(frags, Q):
if len(frags)>0:
    cNodes = []
    zNodes = []
    for node in frags.nodes():
        nodeType, q = node
        if nodeType[0]=='c':
            l = cNodes
        else:
            l = zNodes
        l.append( ( node, frags.node[node]['intensity'] ) )
    nodes = []
    intensities = []
    interactions= []
    costs = []
    for (cz, q), I in cNodes + zNodes:
        intensities.append(I)
        nodes.append( (cz,q) )
        interactions.append(frozenset([(cz,q)]))
        costs.append( Q-1-q )
    for (c, cQ), cI in cNodes:
        for (z, zQ), zI in zNodes:
            if cQ + zQ <= Q-2:
                interactions.append(frozenset([(c,cQ),(z,zQ)]))
                costs.append(Q-1-cQ-zQ)


    A = np.zeros((len(nodes), len(interactions)))
    for i,n in enumerate(nodes):
        for j, interaction in enumerate(interactions):
            if n in interaction:
                A[i,j] = 1
    optim_result = SimplexAlgorithm( costs, A_eq=A, b_eq=intensities, options={"disp": False})
    res = {}
    for inter, I in zip(interactions, optim_result.x):
        if I > 0:
            res[inter] = I
else:
    node = frags.nodes(data=True)[0]
    res = {frozenset(node[0]): node[1]['intensity']}
# return res


fragments[10]
collect_fragments(fragments[2],Q)


for frags in fragments:
    print fragments[frags]
[ collect_fragments(fragments[frags],Q) for frags in fragments]
