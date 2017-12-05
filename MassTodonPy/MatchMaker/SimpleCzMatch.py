from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
import networkx as nx
from networkx import connected_component_subgraphs as connected_components


SimpleNode = namedtuple('N', ['type', 'no', 'bp', 'q'])

# TODO: instantiate all reported variables in the class init
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
        self._make_info()

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
            self.graph.add_edge(C, Z)    

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

    def _optimize(self, G):
        """Match the intensities in a cluster."""
        Q = self.precursor.q
        # lavish: all fragments lose cofragments
        lavish = sum((Q - 1 - N.q) * I for N, I in
                     G.nodes.data('intensity'))
        self.I_lavish += lavish
        if len(G) > 1:  # not trivial
            FG = nx.DiGraph()
            FG.add_node('S')  # start
            FG.add_node('T')  # terminus/sink
            for C, C_intensity in G.nodes.data("intensity"):
                if C.type is 'c':
                    FG.add_node(C)
                    FG.add_edge('S', C, capacity=C_intensity)
                    for Z in G[C]:
                        Z_intensity = G.node[Z]['intensity']
                        FG.add_node(Z)
                        FG.add_edge(C, Z)
                        FG.add_edge(Z, 'T', capacity=Z_intensity)
            total_flow, flows = nx.maximum_flow(FG, 'S', 'T')
            self.I_ETnoD_PTR_fragments += lavish - (Q - 1) * total_flow
            for N in G:
                G.add_edge(N, N)
            # no double count: flows = { start: {end: {value}, .. }, .. }
            for N in flows:
                for M in flows[N]:
                    if N is 'S':  # M is a C fragment
                        G[M][M]['flow'] = G.node[M]['intensity'] - flows[N][M]
                    elif M is 'T':  # N is a Z fragment
                        G[N][N]['flow'] = G.node[N]['intensity'] - flows[N][M]
                    else:  # N is a C and M a Z fragment
                        G[N][M]['flow'] = flows[N][M]
        else:  # trivial case
            self.I_ETnoD_PTR_fragments += lavish
            N, N_intensity = list(G.nodes.data('intensity'))[0]
            G.add_edge(N, N, flow=N_intensity)
            flows = None
        for N, M, intensity in G.edges.data('flow'):
            self.I_ETD_bond[M.bp] += intensity
        return flows

    def _make_info(self):
        """Prepare the information on matches."""
        assert self.I_ETnoD_PTR_fragments >= 0,\
            "The total intensity of ETnoD and PTR on fragments should be non-negative.\
            But it equals {}".format(self.I_ETnoD_PTR_fragments)

        self.I_ETD = sum(v for k, v in self.I_ETD_bond.items())
        self.I_reactions = self.I_ETnoD + self.I_PTR + self.I_ETD
        self.I_unreacted_precursor = self.I_ETnoD_PTR_precursor[0,0]
        self.I_ETD_ETnoD_PTR = self.I_ETnoD + self.I_PTR + self.I_ETD
        # P_ = Probability
        self.P_bond = {k: v/self.I_ETD for k, v in self.I_ETD_bond.items()}
        self.P_reaction = self.I_reactions / (self.I_reactions + self.I_unreacted_precursor)
        self.P_fragmentation = self.I_ETD / self.I_reactions
        self.P_ETD = self.I_ETD / self.I_ETD_ETnoD_PTR
        self.P_ETnoD = self.I_ETnoD / self.I_ETD_ETnoD_PTR
        self.P_PTR = self.I_PTR / self.I_ETD_ETnoD_PTR

    def get_probabilities(self):
        return {k[2:]: v for k, v in self.__dict__.items() if k[0:2] == 'P_'}

    def get_intensities(self):
        return {k[2:]: v for k, v in self.__dict__.items() if k[0:2] == 'I_'}

    def _match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities."""
        self.flows = []
        for G in connected_components(self.graph):
            self.flows.append(self._optimize(G))

