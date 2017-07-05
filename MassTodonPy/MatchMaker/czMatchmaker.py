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

from    collections     import  Counter
import  networkx        as      nx
from    cvxopt          import  matrix, spmatrix, sparse, spdiag, solvers


class czMatchMaker(object):
    '''Virtual class of all matchmakers.'''
    def __init__(   self,
                    MassTodonResults,
                    Q,
                    fasta,
                    accept_nonOptimalDeconv = False,
                    min_acceptEstimIntensity= 100.,
                    verbose = False):
        self.MassTodonResults = MassTodonResults
        self.Q = Q
        self.fasta = fasta
        self.verbose = verbose
        self.min_acceptEstimIntensity = min_acceptEstimIntensity
        self.accept_nonOptimalDeconv = accept_nonOptimalDeconv

    def define_fragment(self, molecule):
        '''Defines what should be considered a node in the c-z matching graphs.'''
        raise NotImplementedError

    def add_edges(self, graph):
        '''Defines what should be considered a node in the c-z matching graphs.'''
        raise NotImplementedError

    def get_graph_analyze_precursor(self):
        '''Generate the graph of pairings, find its connected components, find the number of PTR and ETnoD reactions on precursors.'''
        unreacted_precursors = ETnoDs_on_precursors = PTRs_on_precursors = 0.0
        graph = nx.Graph()
        Q = self.Q
        for res in self.MassTodonResults:
            if (self.accept_nonOptimalDeconv or res['status']=='optimal') and res['status'] != 'ValueError': #TODO what to do otherwise? Nothing for now.
                for mol in res['alphas']:
                    estimate = mol['estimate']
                    if estimate > self.min_acceptEstimIntensity:
                        if mol['molType']=='precursor':
                            q, g = mol['q'], mol['g']
                            if q==Q and g==0:
                                unreacted_precursors = int(estimate)
                            else:
                                ETnoDs_on_precursors += g * estimate
                                PTRs_on_precursors   += (Q-q-g) * estimate
                        else:
                            frag = self.define_fragment(mol)
                            if not frag in graph:
                                graph.add_node( frag, intensity=0 )
                            graph.node[frag]['intensity'] += int(estimate)
        ETnoDs_on_precursors = int(ETnoDs_on_precursors)
        PTRs_on_precursors   = int(PTRs_on_precursors)
        graph = self.add_edges(graph)
        return graph, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors

    def optimize(self):
        '''Perform the optimal pairing of c and z fragments.'''
        raise NotImplementedError

    def get_probs(self, Counts, Probs, tag1, tag2, name1=None, name2=None):
        '''Make two probabilities out of counts and name them properly.'''
        if not name1:
            name1 = tag1
        if not name2:
            name2 = tag2
        total = Counts[tag1]+Counts[tag2]
        if total > 0.0:
            Probs[name1] = float(Counts[tag1])/total
            Probs[name2] = 1.0 - Probs[name1]
        return Probs

    def get_etnod_ptr_probs(self, Counts, Probs):
        Probs = self.get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = self.get_probs(Counts,Probs,'ETnoD','PTR')
        return Probs

    def analyze_counts(self, Counts):
        Counts['total_frags'] = float(sum(Counts[k] for k in Counts if isinstance(k,(int,long))))
        Probs = Counter()
        if Counts['total_frags'] > 0.0:
            for k in Counts:
                if isinstance(k, (int,long)):
                    Probs[k] = float(Counts[k])/Counts['total_frags']
        Counts['total_reactions'] = sum(Counts[k] for k in Counts if k != 'unreacted_precursors')
        Probs = self.get_probs( Counts, Probs, 'unreacted_precursors', 'total_reactions', 'anion_did_not_approach_cation', 'anion_approached_cation' )
        if Counts['total_reactions'] > 0.0:
            Probs['fragmentation'] = float(Counts['total_frags'])/Counts['total_reactions']
            Probs['no fragmentation'] = 1.0 - Probs['fragmentation']
        Counts['ETnoD'] = Counts['ETnoD_frag'] + Counts['ETnoD_precursor']
        Counts['PTR']   = Counts['PTR_frag'] + Counts['PTR_precursor']
        Probs = self.get_etnod_ptr_probs(Counts,Probs)
        return Probs, Counts


    def match(self):
        '''Pair molecules minimizing the number of reactions and calculate the resulting probabilities.'''
        Counts = Counter()

        graph, Counts['ETnoD_precursor'], Counts['PTR_precursor'], Counts['unreacted_precursors'] = self.get_graph_analyze_precursor()

        OptimInfo = []
        for G in nx.connected_component_subgraphs(graph):
            if self.verbose:
                cnt, info = self.optimize(G)
                Counts += cnt
                OptimInfo.append(info)
            else:
                Counts += self.optimize(G)

        Probs, Counts = self.analyze_counts(Counts)

        if self.verbose:
            return Probs, Counts, OptimInfo
        else:
            return Probs, Counts


##################################################################################################




class czMatchMakerBasic(czMatchMaker):
    def define_fragment(self, molecule):
        return molecule['molType'], molecule['q']


    def add_edges(self, graph):
        for C, qC in graph: # adding edges between c and z fragments
            if C[0]=='c':
                for Z, qZ in graph:
                    if Z[0]=='z':
                        bpC = int(C[1:])
                        bpZ = len(self.fasta) - int(Z[1:])
                        if bpC==bpZ and qC + qZ < self.Q-1:
                            graph.add_edge((C,qC),(Z,qZ))
        return graph

    def optimize(self, G):
        '''Finds the minimal number of reactions necessary to explain the MassTodon results.

        Uses the max flow algorithm in all but trivial cases.
        '''
        no_edges_reactions_cnt = sum( (self.Q-1-N[1])*G.node[N]['intensity'] for N in G)
        Counts = Counter()
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
            Counts['reactions'] = no_edges_reactions_cnt - (self.Q-1)*flow_val
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
            Counts['reactions'] = no_edges_reactions_cnt
            G.nodes(data=True)
            N = G.nodes()[0]
            G.add_edge( N, N, flow=G.node[N]['intensity'] )
            flow_val= flows = None

        for N in G:
            for M in G[N]:
                if M[0][0]=='z':
                    fragmented_AA = len(self.fasta) - int(M[0][1:])
                else:
                    fragmented_AA = int(M[0][1:])
                Counts[ fragmented_AA ] += G[N][M]['flow']

        if self.verbose:
            return Counts, (flow_val, flows)
        else:
            return Counts



##################################################################################################



class czMatchMakerIntermediate(czMatchMaker):
    def define_fragment(self, molecule):
        molG = molecule['g']
        molQ = molecule['q']
        if molG == - 1:     # HTR product
            molG += 1
        if molG + molQ == self.Q:# HTR product
            molG -= 1
        return molecule['molType'], molQ, molG

    def etnod_ptr_on_c_z_pairing(self, q0, g0, q1, g1):
        '''Get the number of ETnoD and PTR reactions on a regular edge.'''
        Netnod  = g0 + g1
        Nptr    = self.Q - 1 - g0 - g1 - q0 - q1
        return Netnod, Nptr

    def get_break_point(self, nType ):
        '''Get the amino acid number that was cleft.'''
        if nType[0] == 'c':
            bP = int(nType[1:])
        else:
            bP = len(self.fasta) - int(nType[1:])
        return bP

    def add_edges(self, graph):
        for Ctype, qC, gC in graph: # adding edges between c and z fragments
                if Ctype[0]=='c':
                    for Ztype, qZ, gZ in graph:
                        if Ztype[0]=='z':
                            bpC = self.get_break_point(Ctype)
                            bpZ = self.get_break_point(Ztype)
                            if bpC==bpZ and qC + qZ + gC + gZ <= self.Q-1:
                                ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ )
                                graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD=ETnoD_cnt, PTR=PTR_cnt )
        return graph

    def optimize(self, G):
        if len(G) > 1:
            Jsum = sum( G.node[N]['intensity'] for N in G)
            bP = self.get_break_point( next(G.nodes_iter())[0])
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
            bP  = self.get_break_point(nType)
            TotalPTR   = 0
            TotalETnoD = 0
            TotalFrags = I
            flow_val = flows = None

        Counts = Counter({'ETnoD_frag':TotalETnoD, 'PTR_frag':TotalPTR, bP: TotalFrags})
        if self.verbose:
            return Counts, (flow_val, flows)
        else:
            return Counts

    def get_etnod_ptr_probs(self, Counts, Probs):
        Probs = self.get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = self.get_probs(Counts,Probs,'ETnoD','PTR')
        Probs = self.get_probs(Counts,Probs,'ETnoD_frag','PTR_frag')
        return Probs





##################################################################################################



solvers.options['show_progress'] = False

class czMatchMakerAdvanced(czMatchMakerIntermediate):
    def __init__(self, MassTodonResults, Q, fasta,
                 accept_nonOptimalDeconv  = False,
                 min_acceptEstimIntensity = 100.,
                 verbose = False,
                 L1 = 0.0, L2 = 0.01
        ):
        self.MassTodonResults = MassTodonResults
        self.Q = Q
        self.fasta = fasta
        self.min_acceptEstimIntensity   = min_acceptEstimIntensity
        self.accept_nonOptimalDeconv    = accept_nonOptimalDeconv
        self.L1 = L1
        self.L2 = L2
        self.verbose = verbose


    def add_edges(self, graph):
        for Ctype, qC, gC in graph: # adding edges between c and z fragments
                # add self loops
            N = (Ctype, qC, gC)
            graph.add_edge( N, N, ETnoD_PTR_cnt=self.Q-1-qC)
            if Ctype[0]=='c':
                for Ztype, qZ, gZ in graph:
                    if Ztype[0]=='z':
                        bpC = self.get_break_point(Ctype)
                        bpZ = self.get_break_point(Ztype)
                        if bpC==bpZ and qC + qZ + gC + gZ <= self.Q-1:
                            ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing(qC, gC, qZ, gZ)
                            graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD_PTR_cnt=self.Q-1-qC-qZ, ETnoD=ETnoD_cnt, PTR=PTR_cnt )
        return graph


    def diag(self, val, dim):
        '''Make a sparse identity matrix multiplied by a scalar val.'''
        return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])


    def incidence_matrix(self, G, Jdim, Idim):
        '''Make a sparse incidence matrix of the graph G.'''
        L = spmatrix([], [], [], size=(Jdim,Idim) )
        NodesNo = dict([ (N,i) for i,N in enumerate(G)])
        for j, (N0, N1) in enumerate(G.edges()):
            L[NodesNo[N0],j] = 1
            L[NodesNo[N1],j] = 1
        return L


    def optimize(self, G):
        '''Find regularized max flow.'''
        if len(G) > 1:
            bP= self.get_break_point(next(G.nodes_iter())[0])
            J = matrix([ float(G.node[N]['intensity']) for N in G ])
            Jdim = J.size[0]
            R = matrix([ D['ETnoD_PTR_cnt'] for N0, N1, D in G.edges(data=True) ])
            R = R + self.L1*1.0
            Idim = R.size[0]
            P = self.diag(self.L2*2, Idim)
            L = self.incidence_matrix(G, Jdim, Idim)
            Gmat = self.diag(-1.0, Idim)
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
            bP  = self.get_break_point(nType)
            TotalPTR   = 0
            TotalETnoD = 0
            TotalFrags = I
            S = {'status':'trivial', 'x':I}

        res = Counter({'ETnoD_frag':TotalETnoD, 'PTR_frag':TotalPTR, bP: TotalFrags})
        if self.verbose:
            return res, S
        else:
            return res


######## wrapper


def match_cz_ions(  results_to_pair,
                    Q,
                    fasta,
                    advanced_args,
                    min_acceptEstimIntensity = 0.0,
                    analyzer = 'intermediate',
                    accept_nonOptimalDeconv = True,
                    verbose = False    ):

    if analyzer != 'advanced':
        chosen_analyzer = {
            'basic':            czMatchMakerBasic,
            'intermediate':     czMatchMakerIntermediate,
        }[analyzer](results_to_pair,
                    Q,
                    fasta,
                    accept_nonOptimalDeconv,
                    min_acceptEstimIntensity,
                    verbose     )
    else:
        chosen_analyzer = czMatchMakerAdvanced(
                    results_to_pair,
                    Q,
                    fasta,
                    accept_nonOptimalDeconv,
                    min_acceptEstimIntensity,
                    verbose,
                    **advanced_args )

    return chosen_analyzer.match()
