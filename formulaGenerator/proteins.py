from aminoAcid import AminoAcid
import sys

A = AminoAcid('A')
B = AminoAcid('V')
A += B
A.plot()

class TooShortSequence(Exception): pass

def Protein(fasta):
	if len(fasta) == 0:
		raise TooShortSequence
	P = AminoAcid(fasta[0])
	try:
		for i in xrange(1,len(fasta)):
			P += AminoAcid(fasta[i])
		return P
	except:
		e = sys.exc_info()[0]
		print e
		print fasta[i]


P = Protein('AVTTT'*100)
P.plot()
import igraph as ig
ig.plot(P.G, bbox=(1000,10000), layout=P.G.layout('drl') )


P = Protein('RPKPQQFFGLM')
P = Protein('RP')

A = AminoAcid('F')
A.plot()
A.OH
aas = [ AminoAcid(a) for a in 'RPKPQQFFGLM']

G = aas[0]
G += aas[1]
G += aas[2]
G += aas[3]
G += aas[4]
G += aas[5]
G += aas[6]
G += aas[7]
G += aas[8]
aas[8].plot()
print G
G += aas[9]
print G
[ name for name in G.G.vs['name'] if name[0]=='8' ]

AminoAcid('RPKPQQFFGLM'[8]).plot()

# instead 0 put the index of N_0 or better one of H in NH2.
BFS = A.G.bfsiter(0)
prev= BFS.next()

A.G.vs.find('0_Calpha').index
A.G.vs.find('0_Cbeta').index
A.G.get_eid(7,18)
	
for curr in BFS:
	A.G.get_eid( prev.index, curr.index )
	


print G.OH
G.G.vs.find('8_HO1_2')
G.G.vs.find('8_Hcarboxyl2')
'RPKPQQFFGLM'[9]
aas[9].plot()

aas[8].G['name']

aas[7].OH

A = AminoAcid('Q')
A.plot()
B = AminoAcid('F')
B.plot()
A += B
A.plot()

A.plot()
B.plot()
help(AA.G.vs)


AA.G.vs[AA.H]['name'] = 'Dupa'

print AA.plot()
AA.OH
AA.H
AA.C1
AA.Nalpha

print AA
h = AA.G.vs.find('Calpha').index
AA.G.vs.find('C1').index
AA.G.delete_vertices(h)

AA.plot()
AC = AminoAcid('Y')
AC.plot()

type(AA) == 'aminoAcid.AminoAcid'
type(10)
isinstance( AA, AminoAcid )

print AA

print AA


substanceP = 'RPKPQQFFGLM'
LAA = AA('A')
RAA = T()

AAleft.vs.find('IUPAC'=)

getAttr(RAA,'name')

ig.plot( LAA, **style(G1,label='IUPAC') )
ig.plot( RAA, **style(G2,label='IUPAC') )

from classApproach import Protein

P = Protein('RPKPQQFFGLM')
 P('A')


 def JoinAAs( LAA, RAA ):
	LAA.delOH()
	RAA.delH()
	Nind = LAA.getNind()
	Cind = RAA.getCind()
	return join(LAA, RAA, Nind, Cind)


