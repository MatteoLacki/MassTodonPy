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

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

############################################################

unreacted_precursors = ETnoDs_on_precursors = PTRs_on_precursors = 0.0
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


def etnod_ptr_on_c_z_pairing( q0, g0, q1, g1, Q ):
    '''Get the number of ETnoD and PTR reactions on a regular edge.'''
    Netnod  = g0 + g1
    Nptr    = Q - 1 - g0 - g1 - q0 - q1
    return Netnod, Nptr


def get_break_point( nType, fasta ):
    '''Get the amino acid number that was cleft.'''
    if nType[0] == 'c':
        bP = int(nType[1:])
    else:
        bP = len(fasta) - int(nType[1:])
    return bP


for Ctype, qC, gC in BFG: # adding edges between c and z fragments
    if Ctype[0]=='c':
        for Ztype, qZ, gZ in BFG:
            if Ztype[0]=='z':
                bpC = get_break_point(Ctype, fasta)
                bpZ = get_break_point(Ztype, fasta)
                if bpC==bpZ and qC + qZ + gC + gZ <= Q-1:
                    ETnoD_cnt, PTR_cnt = etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ, Q )
                    BFG.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD=ETnoD_cnt, PTR=PTR_cnt )

fragmentations = set()
for (nT, nQ, nG) in BFG:
    fragmentations.add(get_break_point(nT, fasta))

for (nT, nQ, nG), (mT, mQ, mG) in BFG.edges_iter():
    fragmentations.add(get_break_point(nT, fasta))

eps = + 0.001

fastaLenNoProlines = len(fasta.replace('P',''))
LogProb = dict([ ( i, -log(fastaLenNoProlines) ) for i,f in enumerate(fasta) if f != 'P'])
LogProb['PTR']    = log(.5+eps)
LogProb['ETnoD']  = log(.5-eps)

Graphs = list(nx.connected_component_subgraphs(BFG))

# Counter(map( lambda G: ( len(G),len(G.edges()) ), Graphs ))
# G = [ cc for cc in Graphs if len(cc)==1][2].copy()
# nx.draw_circular(G, with_labels=True, node_size=50 )
# plt.show()


def logBinomial(m,n):
    return lgamma(m+n+1.0)-lgamma(m+1.0)-lgamma(n+1.0)


def etnod_ptr_on_missing_cofragment(nQ, nG, logPetnod, logPptr, Q):
    '''Get the number of ETnoD and PTR reactions on an edges with minimal cost.'''
    if logPetnod > logPptr:
        Netnod  = Q - 1 - nQ
        Nptr    = 0
    else:
        Netnod  = nG
        Nptr    = Q - 1 - nQ - nG
    return Netnod, Nptr


def get_costs(G, Q, LogProb, bP, const=100):
    '''Get the costs.'''
    c  = []
    for C, Z in G.edges_iter():
        if C[0][0]=='z':
            C, Z = Z, C
        (cT, cQ, cG), (zT, zQ, zG) = C, Z
        Netnod      = G.edge[C][Z]['ETnoD']
        Nptr        = G.edge[C][Z]['PTR']
        w_e         = logBinomial(Netnod, Nptr)
        logPptr     = LogProb['PTR']
        logPetnod   = LogProb['ETnoD']
        Cetnod, Cptr = etnod_ptr_on_missing_cofragment(cQ, cG, logPetnod, logPptr, Q)
        Zetnod, Zptr = etnod_ptr_on_missing_cofragment(zQ, zG, logPetnod, logPptr, Q)
        G.node[C]['ETnoD']  = Cetnod
        G.node[C]['PTR']    = Cptr
        G.node[Z]['ETnoD']  = Zetnod
        G.node[Z]['PTR']    = Zptr
        if logPetnod > logPptr:
            W_edge  = (logPptr-logPetnod) * Nptr - (Q-1)*logPetnod + w_e - LogProb[bP]
        else:
            w_cc    = logBinomial(Cetnod, Cptr)
            w_zz    = logBinomial(Zetnod, Zptr)
            W_edge  = -logPptr*(Q-1) + w_e - w_cc - w_zz - LogProb[bP]
        c.append(int(-W_edge * const))
    c = np.array(c)
    return c

def unpaired_cnt(R, G, J):
    return sum( G.node[N][R]*j for N, j in zip(G, J) )

def paired_cnt(R, G, I):
    return sum( (G.edge[N0][N1][R]-G.node[N0][R]-G.node[N1][R])*i for (N0,N1),i in zip(G.edges(),I) )

from scipy.optimize import linprog

def max_weight_flow_simplex(G, Q, LogProb, fasta, verbose=False, const=10000):
    '''Solve one pairing problem.'''
    if len(G) > 1:
        L  = np.array(nx.incidence_matrix(G).todense())
        J  = np.array([G.node[N]['intensity'] for N in G ])
        bP = get_break_point( next(G.nodes_iter())[0], fasta )
        c  = get_costs(G,Q,LogProb,bP)
        I  = linprog(c=c, A_ub=L, b_ub=J )
        TotalFlow = I['x'].sum()
        if LogProb['ETnoD'] < LogProb['PTR']:
            TotalPTR   = unpaired_cnt('PTR', G, J) - (Q-1)*TotalFlow
            TotalETnoD = unpaired_cnt('ETnoD', G, J)
        else:
            pairedPTR  = paired_cnt('PTR', G, I['x'])
            TotalPTR   = unpaired_cnt('PTR', G, J) + pairedPTR
            TotalETnoD = unpaired_cnt('ETnoD', G, J) - (Q-1)*TotalFlow - pairedPTR
        TotalFrags = J.sum() - np.matmul( L, I['x'] ).sum()
        status = I['status']
    else:
        (nType, nQ, nG), Data =  G.nodes(data=True)[0]
        I   = Data['intensity']
        bP  = get_break_point( nType, fasta )
        ETnoD_cnt, PTR_cnt = etnod_ptr_on_missing_cofragment(nQ, nG, LogProb['ETnoD'], LogProb['PTR'], Q)
        TotalFrags = I
        TotalETnoD = I*ETnoD_cnt
        TotalPTR   = I*PTR_cnt
        status = -1
    return Counter({'ETnoD':TotalETnoD, 'PTR':TotalPTR, bP: TotalFrags}), status

%%time
S = Counter()
stati = Counter()
for G in Graphs:
    s, status = max_weight_flow_simplex(G, Q, LogProb, fasta, verbose=False, const=10000)
    S += s
    stati[status] += 1

S['ETnoD'] += ETnoDs_on_precursors
S['PTR'] += PTRs_on_precursors

stati

# Updating Probs
####TODO ADDD THE BLOODY PRECURSORS!!!
LogProb['ETnoD'] = log(S['ETnoD']) - log(S['ETnoD'] + S['PTR'])
LogProb['PTR']   = log(S['PTR'])   - log(S['ETnoD'] + S['PTR'])

TotLogFrag = log(sum( S[s] for s in S if s != 'ETnoD' and s != 'PTR' ))

for s in LogProb:
    if s != 'ETnoD' and s != 'PTR':
        LogProb[s] = log(S[s]) - TotLogFrag

LogProb


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
