# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from    math            import log, lgamma, log10, exp, sqrt
from    collections     import Counter
from    scipy.optimize  import linprog
import  networkx        as nx
import  numpy           as np

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


def get_graphs_analyze_precursors_get_LogProb(MassTodonResults, Q, fasta, eps = 0.0):
    '''Generate the graph of pairings, find its connected components, find the number of PTR and ETnoD reactions on precursors, get the initial logprobabilities of fragments and ETnoD and PTR reactions.'''
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
    ETnoDs_on_precursors    = int(ETnoDs_on_precursors)
    PTRs_on_precursors      = int(PTRs_on_precursors)
    for Ctype, qC, gC in BFG: # adding edges between c and z fragments
        if Ctype[0]=='c':
            for Ztype, qZ, gZ in BFG:
                if Ztype[0]=='z':
                    bpC = get_break_point(Ctype, fasta)
                    bpZ = get_break_point(Ztype, fasta)
                    if bpC==bpZ and qC + qZ + gC + gZ <= Q-1:
                        ETnoD_cnt, PTR_cnt = etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ, Q )
                        BFG.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD=ETnoD_cnt, PTR=PTR_cnt )
    Graphs = list(nx.connected_component_subgraphs(BFG))
    fragmentations = set()
    for (nT, nQ, nG) in BFG:
        fragmentations.add(get_break_point(nT, fasta))
    for (nT, nQ, nG), (mT, mQ, mG) in BFG.edges_iter():
        fragmentations.add(get_break_point(nT, fasta))
    fastaLenNoProlines = len(fasta.replace('P',''))
    LogProb = dict([ ( i, -log(fastaLenNoProlines) ) for i,f in enumerate(fasta) if f != 'P'])
    LogProb['PTR']    = log(.5+eps)
    LogProb['ETnoD']  = log(.5-eps)
    return Graphs, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors, LogProb


def logBinomial(m,n):
    '''The logarithm of the binomial coefficient Bin(m+n, m).'''
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


def get_costs(G, Q, LogProb, bP, const=10000):
    '''Get the costs for the min cost problem.'''
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
    '''Calculate the dot product of estimated MassTodon results and reaction counts for unpaired fragments.'''
    return sum( G.node[N][R]*j for N, j in zip(G, J) )


def paired_cnt(R, G, I):
    '''Calculate the dot product of optimal flow and reaction counts.'''
    return sum( (G.edge[N0][N1][R]-G.node[N0][R]-G.node[N1][R])*i for (N0,N1),i in zip(G.edges(),I) )


def cross_prod_log_binomials_and_J(G, J):
    '''Calculate the sum of logarithms of binomial coefficients for unpaired fragments.'''
    return sum( logBinomial( G.node[N]['ETnoD'], G.node[N]['PTR'] )*j for N, j in zip(G, J) )


# G = [G for G in Graphs if len(G)==1][0]
# G.nodes(data=True)
def max_weight_flow_simplex(G, Q, LogProb, fasta, verbose=False, const=10000):
    '''Solve one pairing problem.'''
    if len(G) > 1:
        L  = np.array(nx.incidence_matrix(G).todense())
        J  = np.array([G.node[N]['intensity'] for N in G ])
        Jsum = J.sum()
        bP = get_break_point( next(G.nodes_iter())[0], fasta )
        c  = get_costs(G,Q,LogProb,bP,const)
        I  = linprog(c=c, A_ub=L, b_ub=J )
        TotalFlow = I['x'].sum()
        unpairedPTR   = unpaired_cnt('PTR', G, J)
        unpairedETnoD = unpaired_cnt('ETnoD', G, J)
        LogLik = I['fun']/float(const) + LogProb[bP]*Jsum + LogProb['ETnoD']*unpairedETnoD
        if LogProb['ETnoD'] < LogProb['PTR']:
            LogLik    += LogProb['PTR']*unpairedPTR + cross_prod_log_binomials_and_J(G, J)
            TotalPTR   = unpairedPTR - (Q-1)*TotalFlow
            TotalETnoD = unpairedETnoD
        else:
            pairedPTR  = paired_cnt('PTR', G, I['x'])
            TotalPTR   = unpairedPTR + pairedPTR
            TotalETnoD = unpairedETnoD - (Q-1)*TotalFlow - pairedPTR
        TotalFrags = Jsum - np.matmul( L, I['x'] ).sum()
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
        LogLik = I * ( LogProb['ETnoD']*ETnoD_cnt + LogProb['PTR']*PTR_cnt + LogProb[bP] )
    return LogLik, Counter({'ETnoD':TotalETnoD, 'PTR':TotalPTR, bP: TotalFrags}), status


def solve_simplex_tasks(Graphs, Q, LogProb, fasta, const=1000):
    '''Solve the linear optimization problems under fixed log probabilities of reactions.'''
    ReactionCount = Counter()
    stati = Counter()
    LogLik= 0.0
    for G in Graphs:
        LogLik_G, s, status = max_weight_flow_simplex(G, Q, LogProb, fasta, verbose=False, const=const)
        LogLik += LogLik_G
        ReactionCount += s
        stati[status] += 1
    return LogLik, ReactionCount, stati


def update_LogProb(LogProb, ReactionCount):
    '''Update the logprobabilities of different reactions.'''
    S = ReactionCount
    LogProb['ETnoD']= log(S['ETnoD']) - log(S['ETnoD'] + S['PTR'])
    LogProb['PTR']  = log(S['PTR'])   - log(S['ETnoD'] + S['PTR'])
    TotLogFrag      = log(sum( S[s] for s in S if s != 'ETnoD' and s != 'PTR' ))
    for s in LogProb:
        if s != 'ETnoD' and s != 'PTR':
            LogProb[s] = log(S[s]) - TotLogFrag
    return LogProb


def L2_distance(PrevLogProb, LogProb):
    return sqrt(sum((LogProb[k]-PrevLogProb[k])**2 for k in LogProb ))


def reaction_analist_advanced(MassTodonResults, Q, fasta, maxIter=100, const=1000, eps = 0.0, crit='logLikDiff', verbose=False, tol=0.01):
    '''Maximize the loglikelihood using coordinate ascent.'''
    Graphs, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors, LogProb = \
        get_graphs_analyze_precursors_get_LogProb(MassTodonResults, Q, fasta, eps)
    i = 0
    if verbose:
        print float(ETnoDs_on_precursors)/(ETnoDs_on_precursors + PTRs_on_precursors)
        print float(PTRs_on_precursors)/(ETnoDs_on_precursors + PTRs_on_precursors)
    LogLikPrev = 0.0
    while True:
        LogLik, ReactionCount, stati = solve_simplex_tasks(Graphs, Q, LogProb, fasta, const)
        ReactionCount['ETnoD'] += ETnoDs_on_precursors
        ReactionCount['PTR']   += PTRs_on_precursors
        PrevLogProb= LogProb
        LogProb    = update_LogProb(LogProb, ReactionCount)
        LogLikDiff = LogLik - LogLikPrev
        LogLikPrev = LogLik
        i += 1
        if verbose:
            print LogLik, stati
            print
        if crit=='logLikDiff':
            stopCond = abs(LogLikDiff)<tol
        if crit=='L2':
            stopCond = L2_distance(PrevLogProb, LogProb)<tol
        if stopCond or i >= maxIter:
            break
    TotalFragmentations = sum(ReactionCount[r] for r in ReactionCount if not r in ('ETnoD', 'PTR') )
    LogProb['no reaction']   = log(unreacted_precursors)-log(unreacted_precursors + TotalFragmentations + ETnoDs_on_precursors + PTRs_on_precursors + ReactionCount['ETnoD'] + ReactionCount['PTR'])
    LogProb['fragmentation'] = log(TotalFragmentations)-log(TotalFragmentations+ReactionCount['ETnoD'] + ReactionCount['PTR'])
    return LogProb
