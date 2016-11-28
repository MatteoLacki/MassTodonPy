import sys
import re
from collections import Counter
try:
  import cPickle as pickle
except:
  import pickle
import igraph as ig


fasta = ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
len(fasta)

modifications = {   ('N', 10) :     Counter({'H': -1, 'O': +2, 'N': +3}),
                    ('Calpha',52) : Counter({'H': -2, 'S': +2, 'N': +2}),
                    ('C',69) :      Counter({'H': -2, 'S': +2, 'N': +2}) }

modifications

G.delete_edges
G.es['Roep']
G.delete_edges(Roep_ne=None)

G.decompose()



for e in G.es(Roep_ne=None):
    print e

for e in G.es:
    print e

position, atomCnt = ('N', 10), Counter({'H': -1, 'O': +2, 'N': +3})
# def checkModification(fasta, position, atomCnt):


def parseModificationKey(tag):
    return {'N':'L', 'Calpha':'C', 'C':'R'}[tag]


def elementContent(G):
    '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
    atomNo = Counter()
    for el in G.vs['elem']:
        atomNo[el] += 1
    return atomNo


if isinstance(fastas, str):
    fastas = [fastas]
assert all( isinstance(f, str) for f in fastas ), 'A fasta you provided was not a string.'

for fasta in fastas:
    for aaNo, aa in enumerate(fasta):


for f in fastas:
    aminos.update(set(f))
assert aminos.issubset(aas), 'One of the fastas contained an unimplemented amino acid.'
aminoAcids = AminoAcids().get()
result = {}

f = fastas[0]
aa = 'F'
# for f in fastas:
lenFasta = len(f)
if lenFasta > 0:
    isSingleton = len(f) == 1
    nTerminus = f[0]
    cTerminus = f[-1]
    f = f[1:-1]
    atomCnt = Counter()
    aminoAcidCounts = Counter(f)
    for aa in aminoAcidCounts:
        G = aminoAcids[aa]['graph']

print(G)
G[]



atomCnt_of_aa = elementContent( G )
for atom in atomCnt_of_aa:
    atomCnt[atom] += atomCnt_of_aa[atom]*aminoAcidCounts[aa]

    # Adding water.
atomCnt['H'] += 2
atomCnt['O'] += 1
atomCnt



modifications

def parseModifications(modifications):


# aminoAcids 			= pickle.load(open('../data/AA.txt', 'rb'))
isoMasses, isoProbs = pickle.load(open('../data/isotopes.txt', 'rb'))
aminoAcidsPattern = re.compile('^(('+')|('.join(list(aminoAcids))+'))+$')
chemFormPattern = re.compile('^((('+')|('.join(list(isoMasses))+'))[0-9]*)+$')
# chemFormPattern = re.compile('^(((O)|(Xe)|(Cs)|(Hg)|(S)|(Ru)|(H)|(Zn)|(Sr)|(Al)|(Sm)|(Zr)|(Ho)|(Ta)|(Pb)|(Te)|(He)|(Ti)|(As)|(Ge)|(Pr)|(U)|(Tl)|(Ir)|(Tm)|(Fe)|(Si)|(Cl)|(Eu)|(Tb)|(W)|(Er)|(P)|(Os)|(K)|(Dy)|(Lu)|(Bi)|(Ga)|(Pt)|(La)|(Be)|(F)|(Yb)|(Kr)|(Cd)|(Mn)|(Ar)|(Cr)|(Se)|(Sb)|(Hf)|(Sc)|(Ca)|(Ba)|(Rb)|(Sn)|(Co)|(Cu)|(Ne)|(Pd)|(In)|(N)|(Au)|(Y)|(Ni)|(Rh)|(C)|(Li)|(Th)|(B)|(Mg)|(Na)|(Pa)|(V)|(Re)|(Nd)|(Br)|(Ce)|(I)|(Ag)|(Gd)|(Nb)|(Mo))[0-9]*)+$')

def isFasta(fasta):
	'''Check if matching a valid fasta sequence'''
	return bool( re.match( aminoAcidsPattern, fasta ) )

def isFormula(formula):
	'''Check if we have a valid chemical formula.'''
	return bool( re.match( chemFormPattern, formula) )

def isPosition(position):
	'''Check if we have a valid position.'''
	return bool( re.match( '^[a-c]\d+$', position) )

def parseFormula(formula):
	'''Return the counter with elemenents as keys and atom counts as numbers.'''
	res = [ (elem, int(count) if len(count)>0 else 1 ) for elem, count in re.findall('([A-Z]+[a-z]*)([0-9]*)', formula) ]
	return Counter(dict(res))

modifications = {}
for arg in sys.argv:
	if '=' in arg:
		preEqual, postEqual = arg.split('=')
		if isPosition(preEqual):
			assert isFormula(postEqual)
			modifications[preEqual] = parseFormula(postEqual)
		elif preEqual == 'fasta':
			assert isFasta(postEqual)
			fasta = postEqual
		else:
			raise IOError

# print(modifications)
# print(fasta)
