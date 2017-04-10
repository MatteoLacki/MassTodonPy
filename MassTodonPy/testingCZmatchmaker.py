%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import  MassTodon
from    Formulator      import  makeFormulas
import  cPickle         as      pickle
import  networkx        as      nx
from    collections     import defaultdict, Counter
from    matplotlib      import collections  as mc
from    cvxopt          import matrix, spmatrix, sparse, spdiag, solvers
from    math            import log, lgamma
import  pylab as pl
import  matplotlib.pyplot as plt
import  numpy           as np
from    graph_compute   import max_cost_flaw

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

############################################################

unreacted_precursors = ETnoDs_on_precursors = PTRs_on_precursors = 0.0
L   = len(fasta)
BFG = nx.Graph()
minimal_estimated_intensity = 100.

for mols, error, status in MassTodonResults:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            if mol['estimate'] > minimal_estimated_intensity: # a work-around the stupidity of the optimization methods
                if mol['molType']=='precursor':
                    if mol['q']==Q and mol['g']==0:
                        unreacted_precursors = mol['estimate']
                    else:
                        ETnoDs_on_precursors+= mol['g'] * mol['estimate']
                        PTRs_on_precursors  += (Q-mol['q']-mol['g']) * mol['estimate']
                else:
                    molG = mol['g']
                    molQ = mol['q']
                    if molG == - 1:     # HTR product
                        molG += 1
                    if molG + molQ == Q:# HTR product
                        molG -= 1
                    frag = (mol['molType'], mol['q'], molG)
                    if not frag in BFG:
                        BFG.add_node( frag, intensity=0 )
                    BFG.node[frag]['intensity'] += int(mol['estimate'])

reactions_on_precursors = ETnoDs_on_precursors + PTRs_on_precursors
# reactions_on_precursors, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors

for C, qC, gC in BFG: # adding edges between c and z fragments
    if C[0]=='c':
        for Z, qZ, gZ in BFG:
            if Z[0]=='z':
                bpC = int(C[1:])
                bpZ = L - int(Z[1:])
                if bpC==bpZ and qC + qZ + gC + gZ <= Q-1:
                    BFG.add_edge( (C,qC,gC), (Z,qZ,gZ))

def etnod_ptr_on_missing_cofragment(nQ, nG, LogProbEtnod, LogProbPtr, Q):
    '''Get the number of ETnoD and PTR reactions on an edges with minimal cost.'''
    if LogProbEtnod > LogProbPtr:
        Netnod  = Q - 1 - nQ
        Nptr    = 0
    else:
        Netnod  = nG
        Nptr    = Q - 1 - nQ - nG
    return Netnod, Nptr

def etnod_ptr_on_c_z_pairing( q0, g0, q1, g1, Q ):
    '''Get the number of ETnoD and PTR reactions on a regular edge.'''
    Netnod  = g0 + g1
    Nptr    = Q - 1 - g0 - g1 - q0 - q1
    return Netnod, Nptr

def get_break_point( nType ):
    '''Get the amino acid number that was cleft.'''
    if nType[0] == 'c':
        bP = int(nType[1:])
    else:
        bP = L - int(nType[1:])
    return bP

fragmentations = set()
for (nT, nQ, nG) in BFG:
    fragmentations.add(get_break_point(nT))

for (nT, nQ, nG), (mT, mQ, mG) in BFG.edges_iter():
    fragmentations.add(get_break_point(nT))


LogProb = dict([ ( i, -log(len(fasta.replace('P',''))) ) for i,f in enumerate(fasta) if f != 'P'])
LogProb['PTR']    = log(.5)
LogProb['ETnoD']  = log(.5)

ccs = list(nx.connected_component_subgraphs(BFG))
# Counter(map( lambda G: ( len(G),len(G.edges()) ), ccs ))

G = [ cc for cc in nx.connected_component_subgraphs(BFG) if len(cc)==8][0]
# nx.draw_circular(G, with_labels=True, node_size=50 )
# plt.show()

def logBinomial(m,n):
    return lgamma(m+n+1.0)-lgamma(m+1.0)-lgamma(n+1.0)


def get_weight(C, Z, LogProb, Q):
    '''Weight for the weighted max flow optimization problem.'''
    (cT, cQ, cG), (zT, zQ, zG) = C, Z
    Netnod, Nptr= etnod_ptr_on_c_z_pairing( cQ, cG, zQ, zG, Q )
    w_e = logBinomial(Netnod, Nptr)

    logPptr     = LogProb['PTR']
    logPetnod   = LogProb['ETnoD']

    bP = get_break_point(cT)
    Cetnod, Cptr = etnod_ptr_on_missing_cofragment(cQ, cG, logPetnod, logPptr, Q)
    Zetnod, Zptr = etnod_ptr_on_missing_cofragment(zQ, zG, logPetnod, logPptr, Q)

    if logPetnod > logPptr:
        W_edge  = (logPptr-logPetnod) * Nptr - (Q-1)*logPetnod + w_e - LogProb[bP]
    else:
        w_cc    = logBinomial(Cetnod, Cptr)
        w_zz    = logBinomial(Zetnod, Zptr)
        W_edge  = -logPptr*(Q-1) + w_e - w_cc - w_zz - LogProb[bP]
    return W_edge


def initialize_flow_graph(G, Q, LogProb):
    '''Construct the flow graph corresponding to one pairing problem.'''

    totalIntensity = sum(G.node[N]['intensity'] for N in G )
    FG = nx.DiGraph()
    FG.add_node('S', demand= -totalIntensity) # start
    FG.add_node('T', demand=  totalIntensity) # terminus/sink
    for C in G:
        if C[0][0]=='c':
            FG.add_node(C)
            FG.add_edge( 'S', C, capacity=G.node[C]['intensity'] )
            for Z in G[C]:
                FG.add_node(Z)
                FG.add_edge( Z, 'T', capacity = G.node[Z]['intensity'] )
                FG.add_edge( C,  Z,   weight  = get_weight(C, Z, LogProb, Q) )
    max_cost_flaw(FG, 'S', 'T', cost="weight", capacity="capacity")
    return FG

# G.nodes(data=True)
# G.edges(data=True)

def initialize_flow_graph_nx(G, Q, LogProb):
    '''Construct the flow graph corresponding to one pairing problem.'''

    totalIntensity = sum(G.node[N]['intensity'] for N in G )
    FG = nx.DiGraph()
    FG.add_node('S', demand= -totalIntensity) # start
    FG.add_node('T', demand=  totalIntensity) # terminus/sink
    FG.add_edge('S','T')
    for C in G:
        if C[0][0]=='c':
            FG.add_node(C)
            FG.add_edge( 'S', C, capacity=G.node[C]['intensity'] )
            for Z in G[C]:
                FG.add_node(Z)
                FG.add_edge(Z, 'T', capacity = G.node[Z]['intensity'])
                weight = int(-get_weight(C,Z,LogProb, Q) * 10000)
                print weight
                FG.add_edge(C, Z, weight = weight)

    print FG.nodes(data=True)
    print FG.edges(data=True)

    minFlaw = nx.min_cost_flow(FG)
    return minFlaw

# nx.draw_circular(G, with_labels=True, node_size=50 )
# plt.show()

from scipy.optimize import linprog

def max_weight_flow_simplex(G, Q, LogProb, verbose=False):
    '''Solve one pairing problem.'''
    L = np.array(nx.incidence_matrix(G).todense())
    J = np.array([G.node[N]['intensity'] for N in G ])
    c = []
    for C, Z in G.edges_iter():
        if C[0][0]=='z':
            C, Z = Z, C
        c.append( int(-get_weight(C,Z,LogProb, Q) * 10000) )
    c = np.array(c)
    simplex_sol = linprog(c=c, A_ub=L, b_ub=J )
    for (N0,N1), flaw in zip(G.edges(), simplex_sol['x']):
        G.edge[N0][N1]['flaw'] = flaw
    if verbose:
        return G, simplex_sol
    else:
        return G



%%time
H = initialize_flow_graph(G.copy(), Q, LogProb)


for C,Z in G.edges_iter():
    if C[0][0]=='z':
        C, Z = Z, C
    print H.edge[C][Z]['flaw']


%%time
H = max_weight_flow_simplex(G, Q, LogProb, verbose=True)

%%time
FGs = [ initialize_flow_graph(G, Q, LogProb) for G in nx.connected_component_subgraphs(BFG)]





%%time
sols = [ max_weight_flow_simplex(G, Q, LogProb, verbose=True) for G in nx.connected_component_subgraphs(BFG)]
Counter( sol['status'] for _,sol in sols)


LogProb


# with open('dupnyGraf.matteo', 'w') as f:
#     pickle.dump(FG, f)
# flowDict = nx.min_cost_flow(FG)


def solve_simplex_step(FGs, Prob, Q):
    '''Finds the maximum a posteriori given probabilities.

    Uses the weighted max flow algorithm in all but trivial cases.
    '''

    max_cost_flaw(FG, 'S', 'T', cost="weight", capacity="capacity")





FG.edges(data=True)


def coordinate_ascent_MLE(FGS, Probs, maxIter=1000):
    for i in xrange(maxIter):
        FGs  = solve_simplex_step(FGs, Prob, Q)
        Prob = solve_analytic_step(FGs, Prob)
    return FGs, Prob





#
#
# fragmentations_no_aas = Counter()
# reactions_on_frags_other_than_fragmentation = 0
# for reactionNo, flow in res:
#     for N in flow:
#         for M in flow[N]:
#             if M=='T':
#                 if N[0][0] == 'c':
#                     breakPoint = int(N[0][1:])
#                 if N[0][0] == 'z':
#                     breakPoint = L-int(N[0][1:])
#                 fragmentations_no_aas[ breakPoint ] += flow[N][M]
#     reactions_on_frags_other_than_fragmentation += reactionNo
#
# ## Testing if all prolines have zero estimates
# # [ fragmentations_no_aas[i] for i,f in enumerate(fasta) if f=='P']
# ## They do.
#
# fragmentations_no_total = sum(fragmentations_no_aas.values())
# fragmentations_no_total
# fragmentations_no_aas
#
# prob_no_reaction = float(no_reactions)/ (no_reactions+reactions_on_precursors+reactions_on_frags_other_than_fragmentation+fragmentations_no_total)
# prob_reaction = 1.0 - prob_no_reaction
#
# prob_fragmentation = float(fragmentations_no_total)/( fragmentations_no_total+reactions_on_frags_other_than_fragmentation+reactions_on_precursors )
#
# prob_no_fragmentation = 1.0 - prob_fragmentation
#
# probs_fragmentation_on_aas = [ float(fragmentations_no_aas[i])/fragmentations_no_total for i in xrange(len(fasta)+1)]
# probs_fragmentation_on_aas
