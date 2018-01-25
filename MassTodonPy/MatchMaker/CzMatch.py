from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
from future.builtins import super

from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch


Node = namedtuple('node', ['type', 'no', 'bp', 'q', 'g'])

class CzMatch(SimpleCzMatch):
    def __init__(self, masstodon):
        """Match c and z ions' intensities.

        Parameters
        ==========
        results : list
            A list of raw results of **MassTodon.run()**.
        precursor : Precursor
            A precursor for the matching problem.

        """
        self._I_ETnoD_fragments = 0
        self._I_PTR_fragments = 0
        super().__init__(masstodon)

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mt, po, cs = molecule._molType_position_cleavageSite()
        return Node(mt, po, cs, molecule.q, molecule.g)
        # mol_type, no, bp, q, g = self._get_node_info(molecule)
        # return Node(mol_type, no, bp, q, g)

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
        Q = self._masstodon.precursor.q
        if C.bp == Z.bp and C.q + Z.q + C.g + Z.g < Q:
            self.graph.add_edge(C, Z, ETnoD= C.g + Z.g,
                                      PTR= Q-1 -C.g -Z.g -C.q -Z.q,
                                      ETnoD_PTR= Q -1 -C.q -Z.q)


## These procedures are false for now. Need to report additionally the min and max intensities of PTR and ETnoD.
## Tests to get PTR_min, PTR_max, ETnoD_min, ETnoD_max failed -> unknown solution.
    # def _get_intensities(self):
    #     """Estimate intensities."""
    #     super()._get_intensities()
    #     self._I_ETnoD_bond = Counter()
    #     self._I_PTR_bond = Counter()
    #     for N, M, intensity in self.graph.edges.data('flow'):
    #         if N is not M:
    #             PTR = self.graph[N][M]['PTR']
    #             self._I_PTR_bond[M.bp] += PTR * intensity
    #             ETnoD = self.graph[N][M]['ETnoD']
    #             self._I_ETnoD_bond[M.bp] += ETnoD * intensity
    #     self._I_PTR_frags = sum(self._I_PTR_bond.values())
    #     self._I_PTR += self._I_PTR_frags
    #     self._I_ETnoD_frags = sum(self._I_ETnoD_bond.values())
    #     self._I_ETnoD += self._I_ETnoD_frags

    # def _get_probabilities(self):
    #     """Estimate probabilities."""
    #     super()._get_probabilities()
    #     if self._I_ETnoD_frags > 0:
    #         self.P_ETnoD_bond = {k: v/self._I_ETnoD_frags for k, v in self._I_ETnoD_bond.items()}
    #     if self._I_PTR_frags > 0:
    #         self.P_PTR_bond = {k: v/self._I_PTR_frags for k, v in self._I_PTR_bond.items()}

    # def _get_branching_ratios(self):
    #     """Estimate branching ratios."""
    #     super()._get_branching_ratios()
    #     if sum(self._I_ETnoD_bond.values()) > 0:
    #         self.branching_ratio['branching_ratio_bond'] = {
    #             self._I_PTR_bond[k]/_I_ETnoD
    #             for k, _I_ETnoD in self._I_ETnoD_bond.items()
    #             if _I_ETnoD > 0}
