from aminoAcid import AminoAcid
import cPickle as pickle

aas = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
aminoAcids = {}
for aa in aas:
    AA = AminoAcid(aa)
    aminoAcids[aa] ={ 'graph' : AA.getGraph(), 'NalphaIDX' : AA.Nalpha(), 'CcarboIDX' : AA.Ccarbo() }

with open('aminoAcids.p','w') as f:
    pickle.dump( aminoAcids, f)
