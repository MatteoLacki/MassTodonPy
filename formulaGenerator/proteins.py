from aminoAcid import AminoAcid

AA = AminoAcid('A')
AA.plot()

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


