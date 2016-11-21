import sys
import re
from collections import Counter
try:
  import cPickle as pickle
except:
  import pickle

aminoAcids 			= pickle.load(open('../data/AA.txt', 'rb'))
isoMasses, isoProbs = pickle.load(open('../data/isotopes.txt', 'rb'))
aminoAcidsPattern 	= re.compile('^(('+')|('.join(list(aminoAcids))+'))+$')
chemFormPattern 	= re.compile('^((('+')|('.join(list(isoMasses))+'))[0-9]*)+$')
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