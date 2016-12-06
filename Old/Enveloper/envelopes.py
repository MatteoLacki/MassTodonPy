from parser import parse
import cPickle as pickle
from collections import Counter

ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

fragments, residues = parse(ubiquitin)
print fragments



AA = pickle.load( open('/Users/matteo/MassSpec/Matteomics/AA.txt', 'rb') )
aas = 'IPP'


atoms = Counter()
for aa in aas:
	atoms += AA[aa]



