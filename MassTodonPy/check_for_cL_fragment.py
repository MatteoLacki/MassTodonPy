%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import  MassTodon
from    Formulator      import  makeFormulas
import  cPickle         as      pickle
import  networkx        as      nx

file_path = '/Users/matteo/Documents/MassTodon/Results/Ubiquitin_ETD_10_ms_1071.matteo'
fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q=8; jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7

with open(file_path, 'rb') as f:
    res = pickle.load(f)

len(fasta)
frags = set()
for re,e,s in res:
    for r in re:
        frags.add(r['molType'])

'z76' in frags
Forms = makeFormulas(fasta=fasta, Q=Q, fragType='cz')
list(Forms.makeMolecules())
frags_made = set()
for t,_,b,_,_ in Forms.makeMolecules():
    frags_made.add(t[0]+str(b))
'z76' in frags_made
res
