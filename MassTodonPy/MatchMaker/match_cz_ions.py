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

from __future__ import absolute_import, division, print_function
from collections import Counter
from future.builtins import super
import networkx as nx
from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers


def make_prob(cnt1, cnt2, ):
    pass


def __get_probs(counts, Probs, tag1, tag2, name1=None, name2=None):
    """Make two probabilities out of counts and name them properly."""
    if not name1:
        name1 = tag1
    if not name2:
        name2 = tag2
    total = counts[tag1]+counts[tag2]
    if total > 0.0:
        Probs[name1] = float(counts[tag1])/total
        Probs[name2] = 1.0 - Probs[name1]
    return Probs


class czMatchMaker(object):
    """Virtual class of all matchmakers.

    Parameters
    ==========
    results : list
        A list of raw results of **MassTodon.run()**.
    precursor : Precursor
        A precursor for the matching problem.
    minimal_intensity : float
        The minimal intenstity of experimental peaks in the deconvolution graph.

    """
    def __init__(self,
                 results,
                 precursor,
                 minimal_intensity=100.0,
                 verbose=False):
        self.results = results
        self.precursor = precursor
        self.verbose = verbose
        self.minimal_intensity = minimal_intensity
        self.intensity = Counter()
        self.get_graph()

    def define_fragment(self, molecule):
        """Defines what should be considered a node in the c-z matching graphs."""
        raise NotImplementedError

    def add_edges_to_graph(self):
        """Defines what should be considered a node in the c-z matching graphs."""
        raise NotImplementedError

    def get_graph(self):
        """Analyze results.

        Returns
        =======
        out : tuple
            The graph of pairings,
            the total intensity of precursors that went through the ETnoD,
            the total intensity of precursors that went through the PTR,
            the total intensity of the unreacted precursors.
        """
        Q = self.precursor.q
        self.graph = nx.Graph()
        for res in self.results:
            if res['status'] is not 'ValueError':
                for mol in res['alphas']:
                    estimate = int(mol['estimate'])
                    if estimate > self.minimal_intensity:
                        mol = mol['molecule']
                        if mol.name is 'precursor':
                            q = mol.q
                            g = mol.g
                            if q is Q and g is 0:
                                self.intensity['unreacted'] = estimate
                            else:
                                self.intensity['ETnoD_precursors'] += g * estimate
                                self.intensity['PTR_precursors'] += (Q - q - g) * estimate
                        else:
                            frag = self.define_fragment(mol)
                            if not frag in self.graph:
                                self.graph.add_node(frag, intensity=0)
                            self.graph.node[frag]['intensity'] += estimate
        self.add_edges_to_graph()

    def optimize(self, G):
        """Perform the optimal pairing of c and z fragments."""
        raise NotImplementedError

    def __get_etnod_ptr_probs(self, Counts, Probs):
        Probs = get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = get_probs(Counts,Probs,'ETnoD','PTR')
        return Probs

    def __analyze_counts(self):
        """Calculate probabilities based on fragment intensities.

        Returns
        =======
        out : tuple
            Probabilities and intensities of reactions.
        """
        self.intensity['total_frags'] = sum(v for k, v in self.intensity.items()
                                            if isinstance(k, (int, long)))
        probs = Counter()
        if self.intensity['total_frags'] > 0.0:
            for k in Counts:
                if isinstance(k, (int,long)):
                    Probs[k] = float(Counts[k])/Counts['total_frags']


        Counts['total_reactions'] = sum(Counts[k] for k in Counts if k != 'unreacted_precursors')
        Probs = get_probs(Counts, Probs, 'unreacted_precursors', 'total_reactions', 'anion_did_not_approach_cation', 'anion_approached_cation' ) # TODO: proportion of untouched precursor
        if Counts['total_reactions'] > 0.0:
            Probs['fragmentation'] = float(Counts['total_frags'])/Counts['total_reactions']
            Probs['no fragmentation'] = 1.0 - Probs['fragmentation']
        Counts['ETnoD'] = Counts['ETnoD_frag'] + Counts['ETnoD_precursor']
        Counts['PTR']   = Counts['PTR_frag'] + Counts['PTR_precursor']
        Probs = self.__get_etnod_ptr_probs(Counts,Probs)
        return Probs, Counts

    def match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities.

        Returns
        =======
        out : tuple
            Probabilities and intensities of reactions.
        """
        OptimInfo = []
        for G in nx.connected_component_subgraphs(self.graph):
            if self.verbose:
                cnt, info = self.optimize(G)
                self.intensity += cnt
                OptimInfo.append(info)
            else:
                self.intensity += self.optimize(G)
        Probs, Counts = self.__analyze_counts(Counts)
        if self.verbose:
            return Probs, Counts, OptimInfo
        else:
            return Probs, Counts


class czMatchMakerBasic(czMatchMaker):
    def define_fragment(self, molecule):
        return molecule.name, molecule.q

    def add_edges_to_graph(self):
        """Add edges between c and z fragments."""
        for C, qC in self.graph:
            if C[0] is 'c':
                for Z, qZ in self.graph:
                    if Z[0] is 'z':
                        bpC = int(C[1:])
                        bpZ = len(self.precursor.fasta) - int(Z[1:])
                        if bpC == bpZ and qC + qZ < self.precursor.q:
                            self.graph.add_edge((C, qC), (Z, qZ))

    def optimize(self, G):
        """Finds the minimal number of reactions necessary to explain the MassTodon results.

        Uses the max flow algorithm in all but trivial cases.
        """
        Q = self.precursor.q
        no_edges_reactions_cnt = sum((Q - 1 - q) * I for (_, q), I in
                                     G.nodes.data('intensity'))
        Counts = Counter()
        if len(G) > 1:  # not trivial
            FG = nx.DiGraph()
            FG.add_node('S') # start
            FG.add_node('T') # terminus/sink
            for C in G:
                if C[0][0] is 'c':
                    Cintensity = G.node[C]['intensity']
                    FG.add_node(C)
                    FG.add_edge('S', C, capacity=Cintensity)
                    for Z in G[C]:
                        Zintensity = G.node[Z]['intensity']
                        FG.add_node(Z)
                        FG.add_edge(C, Z)
                        FG.add_edge(Z, 'T', capacity=Zintensity)
            flow_val, flows = nx.maximum_flow(FG, 'S', 'T')
            Counts['reactions'] = no_edges_reactions_cnt - (Q - 1) * flow_val
            for N in G:
                G.add_edge(N, N)
            for N in flows:
                for M in flows[N]:
                    if N is 'S': # M is a C fragment
                        G.edge[M][M]['flow'] = G.node[M]['intensity']-flows[N][M]
                    elif M is 'T': # N is a Z fragment
                        G.edge[N][N]['flow'] = G.node[N]['intensity']-flows[N][M]
                    else: # N is a C and M a Z fragment
                        G.edge[N][M]['flow'] = flows[N][M]
        else:  # trivial
            Counts['reactions'] = no_edges_reactions_cnt
            G.nodes(data=True)
            N = G.nodes()[0]
            G.add_edge(N, N, flow=G.node[N]['intensity'])
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


class czMatchMakerIntermediate(czMatchMaker):
    def define_fragment(self, molecule):
        molG = molecule.g
        molQ = molecule.q
        mol_name = molecule.name
        if molG == - 1:  # HTR product
            molG += 1
        if molG + molQ == self.precursor.q:  # HTR product
            molG -= 1
        return mol_name, molQ, molG

    def etnod_ptr_on_c_z_pairing(self, q0, g0, q1, g1):
        """Get the number of ETnoD and PTR reactions on a regular edge."""
        Netnod  = g0 + g1
        Nptr    = self.precursor.q - 1 - g0 - g1 - q0 - q1
        return Netnod, Nptr

    def get_break_point(self, nType ):
        """Get the amino acid number that was cleft."""
        if nType[0] == 'c':
            bP = int(nType[1:])
        else:
            bP = len(self.precursor.fasta) - int(nType[1:])
        return bP

    def add_edges_to_graph(self):
        Q = self.precursor.q
        for Ctype, qC, gC in self.graph: # add c-z edges
                if Ctype[0]=='c':
                    for Ztype, qZ, gZ in self.graph:
                        if Ztype[0]=='z':
                            bpC = self.get_break_point(Ctype)
                            bpZ = self.get_break_point(Ztype)
                            if bpC==bpZ and qC + qZ + gC + gZ < Q:
                                ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ )
                                self.graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD=ETnoD_cnt, PTR=PTR_cnt )

    def optimize(self):
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
                    FG.add_edge('S', C, capacity=Cintensity)
                    for Z in G[C]:
                        Zintensity = G.node[Z]['intensity']
                        FG.add_node(Z)
                        FG.add_edge(C, Z)
                        FG.add_edge(Z, 'T', capacity=Zintensity)
            flow_val, flows = nx.maximum_flow(FG,'S','T')
            TotalFrags  = Jsum - flow_val
            TotalETnoD  = 0
            TotalPTR    = 0
            for C, Z in G.edges_iter():
                if C[0][0]=='z':
                    C,Z = Z,C
                # WHAT ABOUT THE UNPAIRED INTENSITIES ???
                #   We cannot tell, how much ETnoD or PTR were there specifically.
                #   All we know, is that there must be a fixed minimal total of both reactions.
                #   So we neglect them in ETnoD and PTR count.
                # BUT WE SHOULD ADD THEM TO OVERALL FRAGMENTATION COUNT, DO WE?
                TotalETnoD  += flows[C][Z]*G.edge[C][Z]['ETnoD']
                TotalPTR    += flows[C][Z]*G.edge[C][Z]['PTR']
        else:
            if self.verbose:
                print('Single node graph.')

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

    def __get_etnod_ptr_probs(self, Counts, Probs):
        Probs = get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = get_probs(Counts,Probs,'ETnoD','PTR')
        Probs = get_probs(Counts,Probs,'ETnoD_frag','PTR_frag')
        return Probs


class czMatchMakerAdvanced(czMatchMakerIntermediate):
    """Advanced intensity matching."""

    def __init__(self,
                 results,
                 precursor,
                 minimal_intensity=100.,
                 verbose=False,
                 L1=0.0,
                 L2=0.01):
        solvers.options['show_progress'] = False
        self.L1 = L1
        self.L2 = L2
        super().__init__(results, precursor, minimal_intensity, verbose)

    def add_edges_to_graph(self):
        for Ctype, qC, gC in self.graph: # adding edges between c and z fragments
                # add self loops
            N = (Ctype, qC, gC)
            self.graph.add_edge( N, N, ETnoD_PTR_cnt=self.Q-1-qC)
            if Ctype[0]=='c':
                for Ztype, qZ, gZ in self.graph:
                    if Ztype[0]=='z':
                        bpC = self.get_break_point(Ctype)
                        bpZ = self.get_break_point(Ztype)
                        if bpC==bpZ and qC + qZ + gC + gZ <= self.Q-1:
                            ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing(qC, gC, qZ, gZ)
                            self.graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD_PTR_cnt=self.Q-1-qC-qZ, ETnoD=ETnoD_cnt, PTR=PTR_cnt)

    def diag(self, val, dim):
        """Make a sparse identity matrix multiplied by a scalar val."""
        return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])

    def incidence_matrix(self, G, Jdim, Idim):
        """Make a sparse incidence matrix of the graph G."""
        L = spmatrix([], [], [], size=(Jdim,Idim) )
        NodesNo = dict([ (N,i) for i,N in enumerate(G)])
        for j, (N0, N1) in enumerate(G.edges()):
            L[NodesNo[N0],j] = 1
            L[NodesNo[N1],j] = 1
        return L

    def optimize(self):
        """Find regularized max flow."""
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


def match_cz_ions(results_to_pair,
                  Q,
                  fasta,
                  advanced_args,
                  minimal_intensity=0.0,
                  analyzer='intermediate',
                  verbose=False):

    if analyzer != 'advanced':
        chosen_analyzer = {'basic': czMatchMakerBasic,
                           'intermediate': czMatchMakerIntermediate
        }[analyzer](results_to_pair,
                    Q,
                    fasta,
                    minimal_intensity,
                    verbose)
    else:
        chosen_analyzer = czMatchMakerAdvanced(results_to_pair,
                                               Q,
                                               fasta,
                                               minimal_intensity,
                                               verbose,
                                               **advanced_args )

    return chosen_analyzer.match()
