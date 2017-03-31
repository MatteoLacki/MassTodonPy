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

from    collections     import defaultdict, Counter
from    matplotlib      import collections  as mc
import  numpy as np
import  pylab as pl
import  matplotlib.pyplot as plt


def update_nators(mol, nominator=0.0, denominator=0.0):
    ETnoDs  = mol['g']
    PTRs    = Q-mol['q']-mol['g']
    nominator   += PTRs*mol['estimate']
    denominator += (ETnoDs+PTRs)*mol['estimate']
    return nominator, denominator

no_reactions = denominator = nominator = 0.0
L = len(fasta)
IDG = nx.Graph() # the intensity division graph
minimal_estimated_intensity = 100.

for mols, error, status in res:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            if mol['estimate'] > minimal_estimated_intensity:
                if mol['molType']=='precursor':
                    if mol['q']==Q and mol['g']==0:
                        no_reactions = mol['estimate']
                    else:
                        nominator, denominator = update_nators(mol, nominator, denominator)
                else:
                    frag = (mol['molType'],mol['q'])
                    IDG.add_node(frag, intensity=mol['estimate'] )
                    IDG.add_edge(frag,frag)

prob_PTR = nominator/denominator
prob_ETnoD = 1.0 - prob_PTR
reactions_on_precursors = denominator

for C, qC in IDG: # adding edges between c and z fragments
    if C[0]=='c':
        for Z, qZ in IDG:
            if Z[0]=='z':
                bpC = int(C[1:])
                bpZ = L - int(Z[1:])
                if bpC==bpZ and qC + qZ < Q-1:
                    IDG.add_edge((C,qC),(Z,qZ))

def min_cost_flow(G, verbose=False):
    '''Finds the minimal number of reactions necessary to explain the MassTodon results.

    Uses the min_cost_flow algorithm.'''
    FG = nx.DiGraph() # flow graph
    FG.add_node('T', demand=0) # we gonna play with integers to eliminate real numbers' imprecisions.
    for N in G:
        intensity = int(G.node[N]['intensity'])
        FG.add_node( N, demand = -intensity )
        FG.node['T']['demand'] += intensity
    for N, M in G.edges_iter():
        if N == M:
            FG.add_edge( N, 'T', weight=Q-1-N[1] )
        else:
            FG.add_node((N,M))
            FG.add_edge( N, (N,M))
            FG.add_edge( M, (N,M))
            FG.add_edge( (N,M), 'T', weight=Q-1-N[1]-M[1])
    res = nx.min_cost_flow_cost(FG)
    if verbose:
        res = (res, nx.min_cost_flow(FG))
    return res

%%time
res = [min_cost_flow(cc, verbose=True) for cc in nx.connected_component_subgraphs(IDG)]

fragmentations_no_aas = Counter()
reactions_on_frags_other_than_fragmentation = 0
for reactionNo, flow in res:
    for N in flow:
        for M in flow[N]:
            if M=='T':
                if N[0][0] == 'c':
                    breakPoint = int(N[0][1:])
                if N[0][0] == 'z':
                    breakPoint = L-int(N[0][1:])
                fragmentations_no_aas[ breakPoint ] += flow[N][M]
    reactions_on_frags_other_than_fragmentation += reactionNo

## Testing if all prolines have zero estimates
# [ fragmentations_no_aas[i] for i,f in enumerate(fasta) if f=='P']
## They do.

fragmentations_no_total = sum(fragmentations_no_aas.values())
fragmentations_no_total
fragmentations_no_aas

prob_no_reaction = float(no_reactions)/ (no_reactions+reactions_on_precursors+reactions_on_frags_other_than_fragmentation+fragmentations_no_total)
prob_reaction = 1.0 - prob_no_reaction

prob_fragmentation = float(fragmentations_no_total)/( fragmentations_no_total+reactions_on_frags_other_than_fragmentation+reactions_on_precursors )

prob_fragmentation
prob_no_fragmentation = 1.0 - prob_fragmentation

probs_fragmentation_on_aas = [ float(fragmentations_no_aas[i])/fragmentations_no_total for i in xrange(len(fasta)+1)]
probs_fragmentation_on_aas
