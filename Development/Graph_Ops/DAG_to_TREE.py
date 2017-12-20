import networkx as nx
import matplotlib.pyplot as plt

DAG = nx.Graph()
DAG.add_node('M0')
options = {'node_color': 'black',
           'node_size': 2,
           'width': 1,
           'with_labels': True,
           'font_size': 20}


for i in range(5):
    I = 'I'+str(i)
    DAG.add_edge('M0', I)
    E = 'E' + str( 0 if i < 2 else 1 )
    DAG.add_edge(I, E)

# nx.draw(DAG, **options); plt.show()

for I in DAG:
    if I[0] is 'I':
        
