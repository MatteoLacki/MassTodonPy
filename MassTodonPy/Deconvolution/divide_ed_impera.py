"""A class to create subproblems for deconvolution.

Should be in principle replaced by a simple C++ code.
Now, it really should."""

import  networkx            as      nx
from    collections         import  defaultdict
import  matplotlib.pyplot   as      plt


class Imperator(object):
    def __init__(self, min_prob = 0.8,
                       isotopic_coverage = .99):
        """Long live the Empire."""
        for prob in (min_prob, isotopic_coverage):
            assert prob > 0.0 and prob <= 1.0
        self.min_prob = min_prob
        self.P        = isotopic_coverage

    def divide(self, molecules, peak_groups):
        """Assign molecules and peak groups to independent problems."""
        G        = nx.Graph()
        min_prob = self.min_prob
        for M in molecules:
            tot_prob = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.P, True):
                E = peak_groups[I_mz]    # index of the group of peaks
                if E > -1:      # -1 denotes emptiness
                    I[(M,E)] += I_prob                 #   E  E  E
                    tot_prob += I_prob                 #    \ | /
            if tot_prob >= min_prob: # planting tree E - M - E in graph G
                G.add_edges_from((M, E, {'prob': P}) for (M, E), P in I.items())
        self.G = G

    def impera(self):
        return nx.connected_component_subgraphs(self.G)

    def plot(self, node_size=.3, show=True):
        """Plot the deconvolution graph."""
        nx.draw(self.G, node_size=node_size)
        if show:
            plt.show()

    def __len__(self):
        return len(self.G)


def divide_ed_impera(molecules,
                     peak_groups,
                     min_prob          = .8,
                     isotopic_coverage = .99):
    imp = Imperator(min_prob, isotopic_coverage)
    imp.divide(molecules, peak_groups)
    return imp