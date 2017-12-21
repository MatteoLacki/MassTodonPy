from collections import Counter, defaultdict
import networkx as nx
import matplotlib.pyplot as plt

from MassTodonPy.Data.Constants import infinity
options = {'node_color': 'black',
           'node_size': 2,
           'width': 1,
           'with_labels': True,
           'font_size': 20}

M = 'M0'

graph = nx.Graph()
graph.add_node(M)
E_2_I = {0: 0, 1: 1, 2: 1, 3:1, 4:2}
for i in range(5):
    I = 'I'+str(i)
    graph.add_node(I, probability=.1, mz=i)
    graph.add_edge('M0', I)
    E = 'E' + str(E_2_I[i])
    graph.add_edge(I, E)
graph.add_edge('I3', 'E2')
graph.add_edge('M0', 'I5')

# nx.draw(graph, **options); plt.show()
"Do not merge isotopologues with no experimental support: their sheer numbers \
indicate that their parent molecule does not belong to the spectrum."

M = 'M0'
visited = {}
I_2_delete = []
for I in graph:
    if I[0] is 'I':
        Gs = frozenset(n for n in graph[I] if n[0] is not 'M')
        if Gs:
            if not Gs in visited:
                visited[Gs] = I
                graph.node[I]['min_mz'] = graph.node[I]['mz']
                graph.node[I]['max_mz'] = graph.node[I]['mz']
                del graph.node[I]['mz']
            else:
                _I = visited[Gs]
                graph.node[_I]['min_mz'] = min(graph.node[_I]['min_mz'],
                                                   graph.node[I]['mz'])
                graph.node[_I]['max_mz'] = max(graph.node[_I]['max_mz'],
                                                   graph.node[I]['mz'])
                graph.node[_I]['probability'] += graph.node[I]['probability']
                I_2_delete.append(I)
graph.remove_nodes_from(I_2_delete)
# nx.draw(graph, **options); plt.show()


def glue_sister_isotopologues(graph):
    """Merge I nodes that share parents.

    The final graph will have I with the id of the smallest I in the group
    of Is that share common E nodes.

    Parameters
    ==========
    graph: networkx.Graph
        The graph consisting of a molecule node M,
        its isotopologues I, and their experimental peaks E.
    """
    visited = {}
    I_2_delete = []
    for I in graph:
        if I[0] is 'I':
            Gs = frozenset(n for n in graph[I] if n[0] is not 'M')
            if Gs:
                if not Gs in visited:
                    visited[Gs] = I
                    graph.node[I]['min_mz'] = graph.node[I]['mz']
                    graph.node[I]['max_mz'] = graph.node[I]['mz']
                    del graph.node[I]['mz']
                else:
                    _I = visited[Gs]
                    graph.node[_I]['min_mz'] = min(graph.node[_I]['min_mz'],
                                                       graph.node[I]['mz'])
                    graph.node[_I]['max_mz'] = max(graph.node[_I]['max_mz'],
                                                       graph.node[I]['mz'])
                    graph.node[_I]['probability'] += graph.node[I]['probability']
                    I_2_delete.append(I)
    graph.remove_nodes_from(I_2_delete)

glue_sister_isotopologues(graph)
# nx.draw(graph, **options); plt.show()
