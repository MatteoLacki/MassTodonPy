"""A class to create subproblems for deconvolution.

Should be in principle replaced by a simple C++ code.
Now, it really should."""

import  networkx            as      nx
from    collections         import  defaultdict, Counter
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
        self.G = nx.Graph()
        for tree in self.divide_iter(molecules, peak_groups):
            self.G.add_edges_from((M_cnt, E, {'prob': P}) for 
                                  (M_cnt, E), P in tree.items())

    def divide_iter(self, molecules, peak_groups):
        for M in molecules:
            tot_prob = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.P, True):
                E = peak_groups[I_mz]    # index of the group of peaks
                if E > -1:      # -1 denotes emptiness
                    I[(M,E)] += I_prob                 #   E  E  E
                    tot_prob += I_prob                 #    \ | /
            if tot_prob >= self.min_prob: # planting tree E - M - E in graph G
                yield I

    def impera_iter(self):
        """Iterate over connected components of the deconvolution graph."""
        return nx.connected_component_subgraphs(self.G)

    def impera(self):
        """List all conected components."""
        self.ccs = list(self.impera_iter())

    def plot(self,
             node_size = 2.0,
             plt_style = 'dark_background',
             show      = True):
        """Plot the deconvolution graph."""
        colors = ['red' if N > 0 else 'blue' for N in self.G]
        plt.style.use(plt_style)
        # red  : peak groups
        # blue : theoretical molecules  
        nx.draw(self.G,
                node_size  = node_size,
                node_color = colors)
        if show:
            plt.show()

    def plot_ccs(self,
                 plt_style = 'dark_background',
                 show      = True):
        """Plot sizes of the small deconvolution graphs."""
        plt.style.use(plt_style)
        edges_cnt = []
        mols_cnt  = []
        peak_groups_cnt = []
        for cc in self.ccs:
            edges_cnt.append(len(cc.edges))
            x = Counter('peaks' if isinstance(N, int) else 'mol' for N in cc)
            mols_cnt.append(x['mol'])
            peak_groups_cnt.append(x['peaks'])
        plt.scatter(mols_cnt, peak_groups_cnt, s = edges_cnt)
        plt.title("Analysis of Connected Components")
        plt.xlabel("Number of Molecules")
        plt.ylabel("Number of Peak Groups")
        if show:
            plt.show()

    def __len__(self):
        return len(self.G)


class ImperatorMagnus(Imperator):
    def divide_iter(self, molecules, peak_groups):
        for M_cnt, M in enumerate(molecules):
            M_cnt = - M_cnt - 1 # it is quicker not to run the __hash__ for M, but use this count
            P_within_groups = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.P, True):
                E = peak_groups[I_mz]
                if E > 0:
                    I[(M_cnt, E)] += I_prob
                    P_within_groups += I_prob
            if P_within_groups >= self.min_prob: 
                # planting tree E - M - E in graph G
                yield I

    def plot_ccs(self,
                 plt_style = 'dark_background',
                 show      = True):
        """Plot sizes of the small deconvolution graphs."""
        plt.style.use(plt_style)
        edges_cnt = []
        mols_cnt  = []
        peak_groups_cnt = []
        for cc in self.ccs:
            edges_cnt.append(len(cc.edges))
            x = Counter(0 if N < 0 else 1 for N in cc)
            mols_cnt.append(x[0])
            peak_groups_cnt.append(x[1])
        plt.scatter(mols_cnt, peak_groups_cnt, s = edges_cnt)
        plt.title("Analysis of Connected Components")
        plt.xlabel("Number of Molecules")
        plt.ylabel("Number of Peak Groups")
        if show:
            plt.show()

    def plot(self,
             node_size = 2.0,
             show      = True):
        """Plot the deconvolution graph."""
        colors = ['red' if N > 0 else 'blue' for N in self.G]
        # red  : peak groups
        # blue : theoretical molecules  
        nx.draw(self.G,
                node_size  = node_size,
                node_color = colors)
        if show:
            plt.show()


# this ultimately selects the default class.
def divide_ed_impera(molecules,
                     peak_groups,
                     min_prob          = .8,
                     isotopic_coverage = .99):
    imp = ImperatorMagnus(min_prob, isotopic_coverage)
    imp.divide(molecules, peak_groups)
    return imp


def draw_connected_component(G, show=True):
    colors= ['red' if N > 0 else 'blue' for N in cc]
    pos   = nx.spring_layout(cc)
    nodes = nx.draw_networkx_nodes(cc,
                                   pos = pos,
                                   node_color = colors)
    edges = nx.draw_networkx_edges(cc,
                                   pos = pos)
    if show:
        plt.show()
