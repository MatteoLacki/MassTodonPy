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


def get_graph_analyze_precursors(MassTodonResults, Q, fasta, minimal_estimated_intensity = 100.):
    '''Generate the graph of pairings, find its connected components, find the number of PTR and ETnoD reactions on precursors.'''
    unreacted_precursors = ETnoDs_on_precursors = PTRs_on_precursors = 0.0
    BFG = nx.Graph()
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
    return BFG, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors


def max_flow(G, fasta):
    if len(G) > 1:
        Jsum = sum( G.node[N]['intensity'] for N in G)
        bP = get_break_point( next(G.nodes_iter())[0], fasta )
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
        TotalFrags  = Jsum - flow_val
        TotalETnoD  = 0
        TotalPTR    = 0
        for C, Z in G.edges_iter():
            if C[0][0]=='z':
                C,Z = Z,C
            TotalETnoD  += flows[C][Z]*G.edge[C][Z]['ETnoD']
            TotalPTR    += flows[C][Z]*G.edge[C][Z]['PTR']
    else:
        (nType, nQ, nG), Data =  G.nodes(data=True)[0]
        I   = Data['intensity']
        bP  = get_break_point( nType, fasta )
        TotalPTR   = 0
        TotalETnoD = 0
        TotalFrags = I
    return Counter({'ETnoD':TotalETnoD, 'PTR':TotalPTR, bP: TotalFrags})


def reaction_analist_intermediate(MassTodonResults, Q, fasta, verbose=False):
    '''Pair molecules minimizing the number of reactions and calculate the resulting probabilities.'''
    Counts = Counter()
    BFG, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors = get_graph_analyze_precursors(MassTodonResults, Q, fasta)
    Counts['ETnoD_precursor']   = ETnoDs_on_precursors
    Counts['PTR_precursor']     = PTRs_on_precursors
    for G in nx.connected_component_subgraphs(BFG):
        Counts += max_flow(G, fasta)
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
        print 'PTR on frags',    Counts['PTR'], 'PTR on prec', Counts['PTR_precursor']
    return Prob
