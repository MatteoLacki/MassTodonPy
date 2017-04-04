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
    MassTodonResults = pickle.load(f)

# from MatchMaker import reaction_analist_basic
# %%time
# reaction_analist_basic(res, fasta, Q) # works perfectly!!!
from    collections     import defaultdict, Counter
from    matplotlib      import collections  as mc
import  pylab as pl
import  matplotlib.pyplot as plt

def update_nators(mol, Q, nominator=0.0, denominator=0.0):
    ETnoDs  = mol['g']
    PTRs    = Q-mol['q']-mol['g']
    nominator   += PTRs*mol['estimate']
    denominator += (ETnoDs+PTRs)*mol['estimate']
    return nominator, denominator

no_reactions = ETnoD_cnt = PTR_cnt = 0.0
L   = len(fasta)
BFG = nx.Graph()
minimal_estimated_intensity = 100.

molTypes = Counter()
for mols, error, status in MassTodonResults:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            molTypes[mol['molType']]+=1

prolines = [i for i,f in enumerate(fasta) if f=='P']
['c'+str(p) in molTypes for p in prolines]
['z'+str(L-p) in molTypes for p in prolines]
# There are no estimates of the proline-cut proteins. Good.

for mols, error, status in MassTodonResults:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            if mol['estimate'] > minimal_estimated_intensity: # a work-around the stupidity of the optimization methods
                if mol['molType']=='precursor':
                    if mol['q']==Q and mol['g']==0:
                        no_reactions = mol['estimate']
                    else:
                        ETnoD_cnt  += mol['g'] * mol['estimate']
                        PTR_cnt    += (Q-mol['q']-mol['g']) * mol['estimate']
                else:
                    frag = ( mol['molType'], mol['q'] )
                    if not frag in BFG:
                        BFG.add_node( frag, intensity=0 ) # intensities will be integers
                    BFG.node[frag]['intensity'] += int(mol['estimate']) # convert the intensities to ints

reactions_on_precursors = ETnoD_cnt + PTR_cnt
prob_PTR   = float(PTR_cnt)/reactions_on_precursors
prob_ETnoD = 1.0 - prob_PTR

for C, qC in BFG: # adding edges between c and z fragments
    if C[0]=='c':
        for Z, qZ in BFG:
            if Z[0]=='z':
                bpC = int(C[1:])
                bpZ = L - int(Z[1:])
                if bpC==bpZ and qC + qZ < Q-1:
                    BFG.add_edge((C,qC),(Z,qZ))

ccs = [ cc for cc in nx.connected_component_subgraphs(BFG) if len(cc)==1]
G = ccs[0]
G.nodes(data=True)
G.edges(data=True)

def minimal_cost(G, Q):
    '''Finds the minimal number of reactions necessary to explain the MassTodon results.

    Uses the max flow algorithm in all but trivial cases.
    '''
        # The number of reactions other than fragmentation
        # that would result from not having any edges between C and Z ions.
    no_edges_reactions_cnt = sum( (Q-1-N[1])*G.node[N]['intensity'] for N in G)

    if len(G)>1:
        FG = nx.DiGraph()
        FG.add_node('S') # start
        FG.add_node('T') # terminus/sink
        for C in G:
            if C[0][0]=='c':
                Cintensity = G.node[C]['intensity']
                FG.add_node(C)
                FG.add_edge( 'S', C, capacity=Cintensity )
                for Z in G[C]:
                    Zintensity = G.node[Z]['intensity']
                    FG.add_node(Z)
                    FG.add_edge(C,Z)
                    FG.add_edge( Z, 'T', capacity=Zintensity )
        flow_val, flows = nx.maximum_flow(FG,'S','T')
        min_cost = no_edges_reactions_cnt - (Q-1)*flow_val
        for N in G:
            G.add_edge(N,N)
        for N in flows:
            for M in flows[N]:
                if N=='S': # M is a C fragment
                    G.edge[M][M]['flow'] = G.node[M]['intensity']-flows[N][M]
                elif M=='T': # N is a Z fragment
                    G.edge[N][N]['flow'] = G.node[N]['intensity']-flows[N][M]
                else: # N is a C and M a Z fragment
                    G.edge[N][M]['flow'] = flows[N][M]
    else:
        min_cost = no_edges_reactions_cnt
        G.nodes(data=True)
        N = G.nodes()[0]
        G.add_edge( N, N, flow=G.node[N]['intensity'] )
    return min_cost, G


ccs = list(nx.connected_component_subgraphs(BFG))
Counter(map(len,ccs))

%%time
res = [ minimal_cost(cc, Q) for cc in nx.connected_component_subgraphs(BFG) if len(cc)>1 ]

# strange_cases = []
# for cc in nx.connected_component_subgraphs(BFG):
#     min_cost, G = minimal_cost(cc, Q)
#     if min_cost < 100:
#         strange_cases.append( (min_cost, G) )
# _,G = strange_cases[0]
# G.nodes(data=True)
# G.edges(data=True)

G = res[0][1]
G.nodes(data=True)
G.edges(data=True)

G[('z12',2)]

fragmentations_no_aas = Counter()
reactions_on_frags_other_than_fragmentation = 0
for reactionNo, G in res:
    reactions_on_frags_other_than_fragmentation += reactionNo
    for N in G:
        for M in G[N]:
            if M[0][0]=='z':
                fragmented_AA = L-int(M[0][1:])
            else:
                fragmented_AA = int(M[0][1:])
            fragmentations_no_aas[ fragmented_AA ] += G[N][M]['flow']

fragmentations_no_total = sum(fragmentations_no_aas.values())

prob_no_reaction = float(no_reactions)/ (no_reactions+reactions_on_precursors+reactions_on_frags_other_than_fragmentation+fragmentations_no_total)
prob_reaction = 1.0 - prob_no_reaction

prob_fragmentation = float(fragmentations_no_total)/( fragmentations_no_total+reactions_on_frags_other_than_fragmentation+reactions_on_precursors )
prob_no_fragmentation = 1.0 - prob_fragmentation

probs_fragmentation_on_aas = [ float(fragmentations_no_aas[i])/fragmentations_no_total for i in xrange(len(fasta)+1)]

results = { 'prob_PTR'          :   prob_PTR,
            'prob_ETnoD'        :   prob_ETnoD,
            'prob_no_reaction'  :   prob_no_reaction,
            'prob_reaction'     :   prob_reaction,
            'prob_fragmentation':   prob_fragmentation,
            'prob_no_fragmentation': prob_no_fragmentation,
            'probs_fragmentation_on_aas':   probs_fragmentation_on_aas }


results 
