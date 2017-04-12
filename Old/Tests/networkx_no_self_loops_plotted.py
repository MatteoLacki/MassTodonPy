import  networkx        as      nx

G = nx.Graph()
G.add_node(0)
G.add_node(1)
G.add_edge(0,0)
G.add_edge(0,1)
G.edges()

pos = nx.spring_layout(G)
nx.draw(G, node_size=10, pos=pos)
nx.draw_networkx_edges(G, pos, width=1.0 )
plt.show()
