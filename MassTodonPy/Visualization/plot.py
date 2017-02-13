import  pylab as pl
from    matplotlib import collections  as mc
import  matplotlib.pyplot as plt
import  networkx as nx
from    collections import defaultdict

def plot_spectrum(spectrum, mz_min=0, mz_max=float('Inf')):
    spec = []
    for i in xrange(spectrum.shape[0]):
        mz, intensity = spectrum[i,0:2]
        if mz >= mz_min and mz <= mz_max:
            spec.append( ( (mz, 0),(mz,intensity) ) )

    lc = mc.LineCollection( spec, linewidths=1 )

    fig, ax = pl.subplots()
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    pl.show()


def plot_connected_component(cc):
    M = []; E = []
    for e in cc.nodes_iter(data=True):
        node_type = e[1]['type']
        if node_type=='M':
            M.append(e[0])
        elif node_type=='E':
            E.append(e[0])

    pos = nx.spring_layout(cc)
    nx.draw_networkx_nodes(cc, pos=pos, nodelist= M, node_color='b', node_size=500, alpha=0.8)
    nx.draw_networkx_nodes(cc, pos=pos, nodelist=E, node_color='r', node_size=20, alpha=0.8)
    nx.draw_networkx_edges(cc, pos, width=1.0 )
    plt.show()


def plot_deconvolution_graph(   G,
                                colors = {'M':'green','I':'blue','G':'red'},
                                sizes  = {'M':40,'I':20,'G':20} ):
    '''Plot the graph for deconvolution problem.

    G       : graph of a deconvolution problem.
    colors  : dictionary mapping tags of graph to colors of nodes
    sizez   : dictionary mapping tags of graph to sizes of nodes
    '''
    pos     = nx.spring_layout(G)
    nodes   = defaultdict(list)
    for N in G:
        nodes[ G.node[N]['type'] ].append(N)
    for Ntype in nodes:
        nx.draw_networkx_nodes(G, pos=pos, nodelist=nodes[Ntype], node_color=colors[Ntype], node_size=sizes[Ntype])
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=.1)
    plt.show()
