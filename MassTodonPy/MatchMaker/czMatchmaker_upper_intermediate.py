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
import  networkx        as     nx
from    collections     import Counter
from    cvxopt      import matrix, spmatrix, sparse, spdiag, solvers
solvers.options['show_progress'] = False


def get_break_point( nType, fasta ):
    '''Get the amino acid number that was cleft.'''
    if nType[0] == 'c':
        bP = int(nType[1:])
    else:
        bP = len(fasta) - int(nType[1:])
    return bP


def etnod_ptr_on_c_z_pairing( q0, g0, q1, g1, Q ):
    '''Get the number of ETnoD and PTR reactions on a regular edge.'''
    Netnod  = g0 + g1
    Nptr    = Q - 1 - g0 - g1 - q0 - q1
    return Netnod, Nptr


def get_graph_analyze_precursors(MassTodonResults, Q, fasta):
    '''Generate the graph of pairings, find its connected components, find the number of PTR and ETnoD reactions on precursors.'''
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
                        N = (mol['molType'], mol['q'], molG)
                        if not N in BFG:
                            BFG.add_node( N, intensity=0 )
                            BFG.add_edge( N, N, ETnoD_PTR_cnt=Q-1-N[1])
                        BFG.node[N]['intensity'] += int(mol['estimate'])
    ETnoDs_on_precursors = int(ETnoDs_on_precursors)
    PTRs_on_precursors   = int(PTRs_on_precursors)
    for Ctype, qC, gC in BFG: # adding edges between c and z fragments
        if Ctype[0]=='c':
            for Ztype, qZ, gZ in BFG:
                if Ztype[0]=='z':
                    bpC = get_break_point(Ctype, fasta)
                    bpZ = get_break_point(Ztype, fasta)
                    if bpC==bpZ and qC + qZ + gC + gZ <= Q-1:
                        ETnoD_cnt, PTR_cnt = etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ, Q )
                        BFG.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD_PTR_cnt=Q-1-qC-qZ, ETnoD=ETnoD_cnt, PTR=PTR_cnt )
    return BFG, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors


def diag(val, dim):
    '''Make a sparse identity matrix multiplied by a scalar val.'''
    return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])


def incidence_matrix(G, Jdim, Idim):
    '''Make a sparse incidence matrix of the graph G.'''
    L = spmatrix([], [], [], size=(Jdim,Idim) )
    NodesNo = dict([ (N,i) for i,N in enumerate(G)])
    for j, (N0, N1) in enumerate(G.edges()):
        L[NodesNo[N0],j] = 1
        L[NodesNo[N1],j] = 1
    return L


def regularized_max_flow(G, fasta, mu=0.0, lam=0.01, verbose=False):
    '''Find regularized max flow.'''
    if len(G) > 1:
        bP= get_break_point( next(G.nodes_iter())[0], fasta )
        J = matrix([ float(G.node[N]['intensity']) for N in G ])
        Jdim = J.size[0]
        R = matrix([ D['ETnoD_PTR_cnt'] for N0, N1, D in G.edges(data=True) ])
        R = R + mu*1.0
        Idim = R.size[0]
        P = diag(lam*2, Idim)
        L = incidence_matrix(G, Jdim, Idim)
        Gmat = diag(-1.0, Idim)
        h = matrix(0.0, size=(Idim, 1))
        S = solvers.qp(P=P, q=R, G=Gmat, h=h, A=L, b=J)
        I = S['x']
        TotalFrags = sum(I)
        TotalPTR   = 0.0
        TotalETnoD = 0.0
        for i, (N0, N1) in enumerate(G.edges_iter()):
            if N0 != N1:
                TotalETnoD  += I[i]*G.edge[N0][N1]['ETnoD']
                TotalPTR    += I[i]*G.edge[N0][N1]['PTR']
    else:
        (nType, nQ, nG), Data = G.nodes(data=True)[0]
        I   = Data['intensity']
        bP  = get_break_point( nType, fasta )
        TotalPTR   = 0
        TotalETnoD = 0
        TotalFrags = I
        S = {'status':'trivial', 'x':I}
    res = Counter({'ETnoD':TotalETnoD, 'PTR':TotalPTR, bP: TotalFrags})
    if verbose:
        return res, S
    else:
        return res


def reaction_analist_upper_intermediate(MassTodonResults, Q, fasta, mu=0.0, lam=0.01, verbose=False):
    '''Pair molecules minimizing the number of reactions and calculate the resulting probabilities.'''
    Counts = Counter()
    BFG, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors = get_graph_analyze_precursors(MassTodonResults, Q, fasta)
    Counts['ETnoD_precursor']   = ETnoDs_on_precursors
    Counts['PTR_precursor']     = PTRs_on_precursors
    stati = []
    for G in nx.connected_component_subgraphs(BFG):
        if verbose:
            res, S = regularized_max_flow(G, fasta, mu, lam, verbose)
            stati.append(S)
        else:
            res = regularized_max_flow(G, fasta, mu, lam, verbose)
        Counts += res
    Prob = Counter()
    TotalReactions = sum(Counts[s] for s in Counts)

    if unreacted_precursors + TotalReactions > 0.0:
        Prob['no reaction'] = float(unreacted_precursors)/(unreacted_precursors + TotalReactions)
        Prob['reaction'] = 1.0 - Prob['no reaction']

    TotalFrags = sum(Counts[s] for s in Counts if isinstance(s, (int,long)) )
    if TotalFrags > 0.0:
        for s in Counts:
            if isinstance(s, (int,long)):
                Prob[s] = float(Counts[s])/TotalFrags

    if TotalReactions > 0.0:
        Prob['fragmentation'] = float(TotalFrags)/TotalReactions
        Prob['no fragmentation'] = 1.0 - Prob['fragmentation']

    TotalETnoD  = Counts['ETnoD']+Counts['ETnoD_precursor']
    TotalPTR    = Counts['PTR']+Counts['PTR_precursor']
    if TotalPTR+TotalETnoD > 0.0:
        Prob['ETnoD'] = float(TotalETnoD)/(TotalPTR+TotalETnoD)
        Prob['PTR']   = 1.0 - Prob['ETnoD']

    if ETnoDs_on_precursors+PTRs_on_precursors > 0.0:
        Prob['ETnoD_prec'] = float(ETnoDs_on_precursors)/(ETnoDs_on_precursors+PTRs_on_precursors)
        Prob['PTR_prec']   = 1.0 - Prob['ETnoD_prec']
        
    if verbose:
        print 'ETnoD on frags',  Counts['ETnoD'], 'ETnoD on prec', Counts['ETnoD_precursor']
        print 'PTR on frags',    Counts['PTR'],   'PTR on prec', Counts['PTR_precursor']
        return Prob, stati
    else:
        return Prob
