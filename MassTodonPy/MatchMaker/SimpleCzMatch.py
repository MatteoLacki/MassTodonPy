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
    minimal_intensity : float
        The minimal intenstity of experimental peaks in the deconvolution graph.

    """
    def __init__(self,
                 results,
                 precursor,
                 minimal_intensity=100.):
        solvers.options['show_progress'] = False
        solvers.options['maxiters'] = 1000
        self.results = results
        self.precursor = precursor
        self.minimal_intensity = minimal_intensity
        # I_ = Intensity
        self.I_ETD_bond = Counter()
        self.I_ETD = 0
        self.I_ETnoD = 0
        self.I_PTR = 0
        self.I_ETnoD_PTR_precursor = Counter()  # len(ETnoD), len(PTR) -> Intensity
        self.I_ETnoD_PTR_fragments = 0
        self.I_lavish = 0
        self._make_graph()
        self._match()
        self._get_intensities()
        self._get_probabilities()
        self._get_branching_ratios()

    def _iter_results(self):
        for res in self.results:
            if res['status'] is not 'ValueError':
                for mol in res['alphas']:
                    estimate = int(mol['estimate'])
                    if estimate > self.minimal_intensity:
                        mol = mol['molecule']
                        yield mol, estimate

    def _get_node_info(self, molecule):
        mol_type = molecule.name[0]
        no = int(molecule.name[1:])
        if mol_type is 'c':
            bp = int(molecule.name[1:])
        else:
            bp = len(self.precursor.fasta) - int(molecule.name[1:])
        return mol_type, no, bp, molecule.q, molecule.g

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mol_type, no, bp, q, g = self._get_node_info(molecule)
        return SimpleNode(mol_type, no, bp, q)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        if C.bp is Z.bp and C.q + Z.q < self.precursor.q:
            self.graph.add_edge(C, Z, ETnoD_PTR=self.precursor.q-1-C.q-Z.q)

    def _add_self_loop(self, N):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        self.graph.add_edge(N, N, ETnoD_PTR=self.precursor.q-1-N.q)

    def _make_graph(self):
        """Prepare the matching graph."""
        Q = self.precursor.q
        self.graph = nx.Graph()
        for mol, estimate in self._iter_results():
            if mol.name is 'precursor':
                self.I_ETnoD += mol.g * estimate
                self.I_PTR += (Q - mol.q - mol.g) * estimate
                self.I_ETnoD_PTR_precursor[mol.g, Q - mol.q - mol.g] = estimate
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
        Q = self.precursor.q
        # lavish: all fragments lose cofragments
        lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
        self.I_lavish += lavish
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
            self.I_ETnoD_PTR_fragments += solution['primal objective']
            for i, (N, M) in enumerate(G.edges()):
                self.graph[N][M]['flow'] = solution['x'][i]
        else:
            self.I_ETnoD_PTR_fragments += lavish
            N, N_intensity = list(G.nodes.data('intensity'))[0]
            self.graph[N][N]['flow'] = N_intensity

    def _get_intensities(self):
        """Estimate intensities."""
        # I_ = Intensity
        assert self.I_ETnoD_PTR_fragments >= 0,\
            "The total intensity of ETnoD and PTR on fragments should be non-negative.\
            But it equals {}".format(self.I_ETnoD_PTR_fragments)
        for N, M, intensity in self.graph.edges.data('flow'):
            self.I_ETD_bond[M.bp] += intensity
        self.I_ETD = sum(v for k, v in self.I_ETD_bond.items())
        self.I_reactions = self.I_ETnoD + self.I_PTR + self.I_ETD
        self.I_unreacted_precursor = self.I_ETnoD_PTR_precursor[0,0]
        self.I_ETD_ETnoD_PTR = self.I_ETnoD + self.I_PTR + self.I_ETD

    def _get_probabilities(self):
        """Estimate probabilities."""
        # P_ = Probability
        if self.I_ETD > 0:
            self.P_ETD_bond = {k: v/self.I_ETD for k, v in self.I_ETD_bond.items()}
        if self.I_reactions + self.I_unreacted_precursor > 0:
            self.P_reaction = self.I_reactions / (self.I_reactions + self.I_unreacted_precursor)
        if self.I_reactions > 0:
            self.P_fragmentation = self.I_ETD / self.I_reactions
        if self.I_ETD_ETnoD_PTR > 0:
            self.P_ETD = self.I_ETD / self.I_ETD_ETnoD_PTR
            self.P_ETnoD = self.I_ETnoD / self.I_ETD_ETnoD_PTR
            self.P_PTR = self.I_PTR / self.I_ETD_ETnoD_PTR

    def _get_branching_ratios(self):
        """Estimate branching ratios."""
        self.branching_ratio = {}
        if self.I_ETnoD and self.I_PTR:
            self.branching_ratio['branching_ratio'] = self.I_PTR / self.I_ETnoD

    def get_probabilities(self):
        return {k[2:]: v for k, v in self.__dict__.items() if k[0:2] == 'P_'}

    def get_intensities(self):
        return {k[2:]: v for k, v in self.__dict__.items() if k[0:2] == 'I_'}

    def get_branching_ratios(self):
        return self.branching_ratio

    def _match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities."""
        for G in connected_components(self.graph):
            self._optimize(G)
