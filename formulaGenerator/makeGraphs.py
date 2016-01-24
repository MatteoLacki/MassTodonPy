from aminoAcid import AminoAcid
import cPickle as pickle

aas 	= ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AAS 	= {}
for aa in aas:
	A = AminoAcid(aa)
	AAS[aa] = ( A.getLeft(), A.getRight(), A.getGraph() )

with open('aminoAcidGraphs.txt','w') as f:
	pickle.dump( AAS, f)