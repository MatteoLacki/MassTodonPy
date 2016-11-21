from IsoSpecPy import IsoSpecPy
from math import exp
from formulaGenerator.fragments import ubiquitin, fasta2atomCount

	# Getting atomic composition of a fasta.
ubi_atomCnt 	= fasta2atomCount(ubiquitin)
ubi_atomCnt_str = "".join( a+str(ubi_atomCnt[a]) for a in ubi_atomCnt ) 

import timeit

def f():
	i = IsoSpecPy.IsoSpec.IsoFromFormula( ubi_atomCnt_str, 0.5)

N 	= 100
time= timeit.repeat('f()', 'from __main__ import f', number=1,repeat=N)
max(time)*20000/60
