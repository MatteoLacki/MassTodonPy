These have to be generated just once.
Counter(map( lambda G: ( len(G),len(G.edges()) ), Graphs ))
G = [ cc for cc in Graphs if len(cc)==1][2].copy()
nx.draw_circular(G, with_labels=True, node_size=50 )
plt.show()
