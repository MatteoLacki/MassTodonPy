from aminoAcid import AminoAcid
try:
  import cPickle as pickle
except:
  import pickle

import sys


def generateAminoAcids():
	'''Based on graph structure generate counters.'''
	if sys.version_info > (3, 0):
		print('Use only Python < 3.0')
	else:
		aas = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
		aminoAcids = {}
		for aa in aas:
		    AA = AminoAcid(aa)
		    aminoAcids[aa] ={ 'graph' : AA.getGraph(), 'NalphaIDX' : AA.Nalpha(), 'CcarboIDX' : AA.Ccarbo() }
		with open('aminoAcids.p','w') as f:
		    pickle.dump( aminoAcids, f)

with open('aminoAcids.p', 'rb') as f:
    d = pickle.load(f)
