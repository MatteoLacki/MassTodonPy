from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers
import networkx as nx
from networkx import connected_component_subgraphs as connected_components

from MassTodonPy.Data.Constants import eps, infinity

SimpleNode = namedtuple('simple_node', ['type', 'no', 'bp', 'q'])

def diag(val, dim):
    """Make a sparse identity matrix multiplied by a scalar val."""
    return spdiag([spmatrix(val,[0],[0]) for i in range(dim)])

def incidence_matrix(graph, row_cnt, col_cnt):
    """Make a sparse incidence matrix of the graph G."""
    L = spmatrix([], [], [], size=(row_cnt, col_cnt) )
    NodesNo = dict([ (N,i) for i,N in enumerate(graph)])
    for j, (N0, N1) in enumerate(graph.edges()):
        L[NodesNo[N0],j] = 1
        L[NodesNo[N1],j] = 1
    return L


class SimpleCzMatch(object):
    """Match c and z ions' intensities neglecting the quenched charge.

    Parameters
    ==========
    results : list
        A list of raw results of **MassTodon.run()**.
    precursor : Precursor
        A precursor for the matching problem.

    """
    def __init__(self,
                 results,
                 precursor):
        solvers.options['show_progress'] = False
        solvers.options['maxiters'] = 1000
        self._results = results
        self._precursor = precursor
        # _I_ = Intensity
        self._I_ETD_bond = Counter()
        self._I_ETD = 0
        self._I_ETnoD = 0
        self._I_PTR = 0
        self._I_ETnoD_PTR_precursor = Counter()  # len(ETnoD), len(PTR) -> Intensity
        self._I_ETnoD_PTR_fragments = 0
        self._I_lavish = 0
        self._make_graph()
        self._match()
        self.intensities = self._get_intensities()
        self.probabilities = self._get_probabilities()
        self.branching_ratio = self._I_PTR / self._I_ETnoD if self._I_ETnoD else None

    def _get_node_info(self, molecule):
        mol_type = molecule.name[0]
        no = int(molecule.name[1:])
        if mol_type is 'c':
            bp = int(molecule.name[1:])
        else:
            bp = len(self._precursor.fasta) - int(molecule.name[1:])
        return mol_type, no, bp, molecule.q, molecule.g

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mol_type, no, bp, q, g = self._get_node_info(molecule)
        return SimpleNode(mol_type, no, bp, q)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        if C.bp is Z.bp and C.q + Z.q < self._precursor.q:
            self.graph.add_edge(C, Z, ETnoD_PTR=self._precursor.q-1-C.q-Z.q)

    def _add_self_loop(self, N):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        self.graph.add_edge(N, N, ETnoD_PTR=self._precursor.q-1-N.q)

    def _make_graph(self):
        """Prepare the matching graph."""
        Q = self._precursor.q
        self.graph = nx.Graph()
        for mol, estimate in self._results:
            if mol.name is 'precursor':
                self._I_ETnoD += mol.g * estimate
                self._I_PTR += (Q - mol.q - mol.g) * estimate
                self._I_ETnoD_PTR_precursor[mol.g, Q - mol.q - mol.g] = estimate
            else:
                frag = self._get_node(mol)
                if not frag in self.graph:
                    self.graph.add_node(frag, intensity=0)
                    self.graph.node[frag]['intensity'] += estimate
        for C in self.graph:
            if C.type is 'c':
                for Z in self.graph:
                    if Z.type is 'z':
                        self._add_edge(C, Z)
        for N in self.graph:
            self._add_self_loop(N)

    def _optimize(self, G):
        """Match the intensities in a cluster.

        Decorates the self.graph with flow values.
        """
        Q = self._precursor.q
        # lavish: all fragments lose cofragments
        lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
        self._I_lavish += lavish
        if len(G) > 1:
            intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
            costs = matrix([float(ETnoD_PTR) for N,M, ETnoD_PTR 
                            in G.edges.data('ETnoD_PTR')])
            edges_cnt = G.size()  # number of c-z pairings
            equalities = incidence_matrix(G, len(intensities), edges_cnt)
            inequalities = diag(-1.0, edges_cnt)
            upper_bounds = matrix([0.0] * edges_cnt)
            primalstart = {}
            primalstart['x'] = matrix([0.0] * edges_cnt)
            primalstart['s'] = matrix([eps] * len(upper_bounds))
            solution = solvers.conelp(c=costs,
                                      G=inequalities,
                                      h=upper_bounds,
                                      A=equalities,
                                      b=intensities,
                                      primalstart=primalstart)
            self._I_ETnoD_PTR_fragments += solution['primal objective']
            for i, (N, M) in enumerate(G.edges()):
                self.graph[N][M]['flow'] = solution['x'][i]
        else:
            self._I_ETnoD_PTR_fragments += lavish
            N, N_intensity = list(G.nodes.data('intensity'))[0]
            self.graph[N][N]['flow'] = N_intensity

    def _get_intensities(self):
        """Estimate intensities."""
        # _I_ = Intensity
        assert self._I_ETnoD_PTR_fragments >= 0,\
            "The total intensity of ETnoD and PTR on fragments should be non-negative.\
            But it equals {}".format(self._I_ETnoD_PTR_fragments)

        for N, M, intensity in self.graph.edges.data('flow'):
            self._I_ETD_bond[M.bp] += intensity

        self._I_ETD = sum(v for k, v in self._I_ETD_bond.items())
        self._I_reactions = self._I_ETnoD + self._I_PTR + self._I_ETD
        self._I_unreacted_precursor = self._I_ETnoD_PTR_precursor[0,0]
        self._I_ETD_ETnoD_PTR = self._I_ETnoD + self._I_PTR + self._I_ETD
        return {k[3:]: v for k, v in self.__dict__.items() if k[0:3] == '_I_'}

    def _get_probabilities(self):
        """Estimate probabilities."""
        # _P_ = Probability
        if self._I_ETD > 0:
            self._P_ETD_bond = {k: v/self._I_ETD for k, v in self._I_ETD_bond.items()}

        if self._I_reactions + self._I_unreacted_precursor > 0:
            self._P_reaction = self._I_reactions / (self._I_reactions + self._I_unreacted_precursor)

        if self._I_reactions > 0:
            self._P_fragmentation = self._I_ETD / self._I_reactions

        if self._I_ETD_ETnoD_PTR > 0:
            self._P_ETD = self._I_ETD / self._I_ETD_ETnoD_PTR
            self._P_ETnoD = self._I_ETnoD / self._I_ETD_ETnoD_PTR
            self._P_PTR = self._I_PTR / self._I_ETD_ETnoD_PTR
        return {k[3:]: v for k, v in self.__dict__.items() if k[0:3] == '_P_'}

    def _match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities."""
        for G in connected_components(self.graph):
            self._optimize(G)
