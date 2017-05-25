from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

mol = substanceP.copy()

cutOff = 100; topPercent = .999; max_times_solve=30
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = True

M = MassTodon(  fasta           = mol['fasta'],
                precursorCharge = mol['Q'],
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = mol['modifications']  )

M.readSpectrum( spectrum        = mol['spectrum'],
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent  )

M.prepare_problems(M_minProb)

Results = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )

SFG = Results[0]['SFG']
GS = [G for G in SFG.nodes(data=True) if G[1]['type']=='G']
A = GS[0][1]['mz']
A.begin
A.end


# SFG = next(M.problems)
# SFG = next(M.problems)
# SFG.edges(data=True)
# SFG.nodes(data=True)
SFG.node['I28']['mz']
# Results = M.run(solver = solver, method  = method,
#                 max_times_solve = max_times_solve,
#                 L1_x = L1_x,        L2_x=L2_x,
#                 L1_alpha=L1_alpha,  L2_alpha=L2_alpha,
#                 verbose = verbose )
#
inf = 10e100

from    intervaltree    import Interval as II, IntervalTree as Itree
from    interval        import interval as III, inf

x = II(3.5, 3.5)
x.begin, x.end

K = Itree((II(3.4, 5., 'a'),  II(3.8, inf, 'b') ))
K = Itree((II(3.4, 5.),  II(3.8, inf), II(5.5, 9) ))
L = Itree()
L.addi(3.5, 4.5)

K
L

K - L
II(1,2).begin

K[ II(3.5,3.9) ]
K[4]

for M in K[3.9]:
    print M.data


x = Counter({'a':1,'b':2})
for k in x.items():
    print k


III[3,4] & III[3.4,inf]

K = Itree((II(3.4, 5.),  II(3.8, inf) ))
K.split_overlaps()
K.addi(3,5)
x = defaultdict(lambda: inf)
x['k'] = min(x['k'],3)
x['k']



def add_zero_intensity_G_nodes(P):
    '''Pair unpaired isotopologue peaks I with zero-intensity experimental groups G.
    '''
    cnts = Counter(P.node[N]['type'] for N in P)
    Gcnt = cnts['G']+1
    newGnodes = []
    newGIedges= []
    for I in P:
        if P.node[I]['type'] == 'I':
            if len(P[I]) == 1:
                G = 'G' + str(Gcnt)
                newGnodes.append( (G,{'intensity':0.0, 'type':'G'}) )
                newGIedges.append( (G,I) )
                Gcnt += 1
    P.add_nodes_from(newGnodes)
    P.add_edges_from(newGIedges)

nx.draw(SFG, with_labels = True)
plt.show()

def create_G_nodes(SFG):
    '''Collect experimental peaks into groups of experimental data G.

    This is done without any loss in resolution.'''
    E2remove = []
    Gs = Counter()
    for E in SFG:
        if SFG.node[E]['type'] == 'E':
            G = frozenset(SFG[E])
            Gs[G] += SFG.node[E]['intensity']
            E2remove.append(E)
    SFG.remove_nodes_from(E2remove)
    for Gcnt, (Is, G_intensity) in enumerate(Gs.items()):
        G = 'G' + str(Gcnt)
        SFG.add_node(G, intensity=G_intensity, type='G')
        for I in Is:
            SFG.add_edge(I, G)
    add_zero_intensity_G_nodes(SFG)



nx.draw(SFG, with_labels = True)
plt.show()
# create_G_nodes(SFG)
# nx.draw(SFG, with_labels = True)
# plt.show()
