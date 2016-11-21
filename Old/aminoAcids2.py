# AAS = {}
# for aa in aas:
#     AAS[aa] = AminoAcid(aa)

# print(AAS)
# print(AAS['A'].G)
# ig.plot(AAS['A'].G)
# print(AminoAcid('V').getGraph())


    # def __str__(self):
    #   return self.G.__str__()

    # def __repr__(self):
    #   return self.G.__repr__()

# def addBackbone(A):
#   B = Backbone()
#   C = join( A, B, vertex_attributes=['elem','name'], edge_attributes=['Roep'] )
#   C += ('Calpha', 'Cbeta')
#   return C

# def makeAtoms(atomCnt):
#   vAtt = {'name': [], 'elem': [] }
#   for atom in atomCnt:
#       vAtt['elem'].extend( [ atom ]*atomCnt[atom] )
#       for aC in range(atomCnt[atom]):
#           vAtt['name'].append( atom+str(aC) )
#       G = ig.Graph( sum( atomCnt[a] for a in atomCnt ), vertex_attrs=vAtt )
#   return G

# def Backbone():
#   B = makeAtoms({ 'C':2, 'H':2, 'N':1, 'O':1 })
#   B.add_edge( 'N0', 'C0', Roep='cz' )
#   B.add_edge( 'C0', 'C1', Roep='ax' )
#   B = B + ('H0','N0')+('C0','H1')+('C1','O0')+('C1','O0')
#   B.vs['name'] = ['HNalpha', 'Halpha', 'Calpha', 'Ccarbo', 'Ocarbo', 'Nalpha']
#   return B

# def plot(G,target=''):
#   elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
#   visual_style = {}
#   visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
#   visual_style['vertex_label']= G.vs['name']
#   visual_style['edge_label']  = G.es['Roep']
#   if len(target) == 0:
#       ig.plot( G, **visual_style )
#   else:
#       ig.plot( G, target=target, **visual_style )

# Alternative tryptophan.
# W = self.makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
# ('C0','H0')+('C0','H1')+('C0','C1')+\
# ('C1','C2')+('C1','C3')+('C1','C3')+\
# ('C2','C4')+('C2','C5')+('C2','C5')+\
# ('C3','N0')+('C3','H2')+\
# ('N0','H3')+('N0','C5')+\
# ('C5','C6')+\
# ('C6','H4')+('C6','C7')+('C6','C7')+\
# ('C4','C8')+('C4','C8')+('C4','H5')+\
# ('C7','H6')+('C7','C8')+\
# ('C8','H7')
