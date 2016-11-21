try:
  import cPickle as pickle
except:
  import pickle

from collections import Counter
from proteinAlternative import roepstorffy, makeFragments, ubiquitin

fasta = ubiquitin
fasta = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"

def mF(fasta):
    fragments = []
    for f in makeFragments( fasta,['cz'] ):
        f[-1].update(Counter({'C':0,'H':0,'N':0,'O':0,'S':0}))
        f.append( roepstorffy(f,fasta) )
        fragments.append(f)
    return fragments

pickle.dump( mF(fasta), open( "test.p", "wb" ) )

fragments = pickle.load(open('/Users/matteo/Documents/Science/MassTodon/MassTodonPy/formulaGenerator/test.p', 'rb'))


### Comparing Residuals # Perfect match
# aminoAcids  = pickle.load(open('aminoAcids.p', 'rb'))
# def elementContent(G):
#     '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
#     atomNo = Counter()  
#     for el in G.vs['elem']:
#         atomNo[el] += 1
#     return atomNo

# aminoAcids = dict([ (aa, elementContent(aminoAcids[aa]['graph'])) for aa in aminoAcids])
