from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
from future.builtins import super

from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch


Node = namedtuple('node', ['type', 'no', 'bp', 'q', 'g'])

class CzMatch(SimpleCzMatch):
    def __init__(self,
                 results,
                 precursor,
                 minimal_intensity=100.):
        """Match c and z ions' intensities.

        Parameters
        ==========
        results : list
            A list of raw results of **MassTodon.run()**.
        precursor : Precursor
            A precursor for the matching problem.
        minimal_intensity : float
            The minimal intenstity of experimental peaks in the deconvolution graph.

        """
        self.I_ETnoD_fragments = 0
        self.I_PTR_fragments = 0
        super().__init__(results, precursor, minimal_intensity)

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mol_type, no, bp, q, g = self._get_node_info(molecule)
        return Node(mol_type, no, bp, q, g)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        # N_PTR = precursor.q - 1 - C.q - Z.q - C.g - Z.g
        #   N_PTR >= 0
        #       N_PTR = precursor.q - 1 - C.q - Z.q - C.g - Z.g
        #       precursor.q - 1 - C.q - Z.q - C.g - Z.g >= 0
        #       precursor.q > C.q + Z.q + C.g + Z.g
        #   N_PTR <= precursor.q - 1
        #       C.q + Z.q + C.g + Z.g >= 0  # automatically
        # N_ETnoD = C.g + Z.g
        #   N_ETnoD >= 0
        #       C.g + Z.g >= 0              # automatically
        #   N_ETnoD < precursor.q
        #       C.g + Z.g < precursor.q,
        #       but
        #       C.q + Z.q + C.g + Z.g < precursor.q
        if C.bp == Z.bp and C.q + Z.q + C.g + Z.g < self.precursor.q:
            self.graph.add_edge(C, Z, ETnoD=C.g + Z.g, 
                                      PTR=self.precursor.q - 1 - C.g - Z.g - C.q - Z.q,
                                      ETnoD_PTR=self.precursor.q - 1 - C.q - Z.q)

    def _get_intensities(self):
        """Estimate intensities."""
        super()._get_intensities()
        self.I_ETnoD_bond = Counter()
        self.I_PTR_bond = Counter()
        for N, M, intensity in self.graph.edges.data('flow'):
            if N is not M:
                PTR = self.graph[N][M]['PTR']
                self.I_PTR_bond[M.bp] += PTR * intensity
                ETnoD = self.graph[N][M]['ETnoD']
                self.I_ETnoD_bond[M.bp] += ETnoD * intensity
        self.I_PTR_frags = sum(self.I_PTR_bond.values())
        self.I_PTR += self.I_PTR_frags
        self.I_ETnoD_frags = sum(self.I_ETnoD_bond.values())
        self.I_ETnoD += self.I_ETnoD_frags

    def _get_probabilities(self):
        """Estimate probabilities."""
        super()._get_probabilities()
        if self.I_ETnoD_frags > 0:
            self.P_ETnoD_bond = {k: v/self.I_ETnoD_frags for k, v in self.I_ETnoD_bond.items()}
        if self.I_PTR_frags > 0:
            self.P_PTR_bond = {k: v/self.I_PTR_frags for k, v in self.I_PTR_bond.items()}
    
    def _get_branching_ratios(self):
        """Estimate branching ratios."""
        super()._get_branching_ratios()
        if sum(self.I_ETnoD_bond.values()) > 0:
            self.branching_ratio['branching_ratio_bond'] = {
                self.I_PTR_bond[k]/I_ETnoD 
                for k, I_ETnoD in self.I_ETnoD_bond.items()
                if I_ETnoD > 0}
