from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
from future.builtins import super

from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch


Node = namedtuple('N', ['type', 'no', 'bp', 'q', 'g'])

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
        super.__init__(results, precursor, minimal_intensity)

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mol_type, no, bp, q, g = self._get_node_info(molecule)
        return Node(mol_type, no, bp, q, g)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        if C.bp is Z.bp and C.q + Z.q + C.g + Z.g < self.precursor.q:
            ETnoD = C.g + Z.g
            PTR = self.precursor.q - 1 - C.g - Z.g - C.q - Z.q
            self.graph.add_edge(C, Z, ETnoD=ETnoD, PTR=PTR)

    def _optimize(self, G):
        if len(G) > 1:
            passycon            

