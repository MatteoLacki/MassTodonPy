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
# modifications = {   ('N',2) :       {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': 1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': 2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': 2, 'S': +2, 'N': +200} }
# Forms = makeFormulas(fasta=fasta, Q=Q, fragType='cz')
# M = MassTodon(  fasta           = fasta,
#                 precursorCharge = Q,
#                 precDigits      = precDigits,
#                 jointProbability= jP,
#                 mzPrec          = mzPrec )
# path  = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
# path += 'Ubiquitin_ETD_10 ms_1071.mzXML'
# cutOff = 100; topPercent = .999
# M.readSpectrum(path=path, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
# M.prepare_problems(M_minProb)
# mu=1e-5; lam=0.0; nu=0.001
# %%time
# res = M.run(solver='sequential', method='MSE', mu=mu, lam=lam, nu=0.001)
# with open(file_path, 'w') as f:
    # pickle.dump(res, f)
with open(file_path, 'rb') as f:
    MassTodonResults = pickle.load(f)

from    collections     import defaultdict, Counter
from    matplotlib      import collections  as mc
import  pylab as pl
import  matplotlib.pyplot as plt

no_reactions = ETnoD_cnt = PTR_cnt = 0.0
L   = len(fasta)
BFG = nx.Graph()
minimal_estimated_intensity = 100.

for mols, error, status in MassTodonResults:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            if mol['estimate'] > minimal_estimated_intensity: # a work-around the stupidity of the optimization methods
                if mol['molType']=='precursor':
                    if mol['q']==Q and mol['g']==0:
                        no_reactions = mol['estimate']
                    else:
                        ETnoD_cnt  += mol['g'] * mol['estimate']
                        PTR_cnt    += (Q-mol['q']-mol['g']) * mol['estimate']
                else:
                    molG = max(mol['g'],0) # glue together HTR and ETD
                    frag = (mol['molType'], mol['q'], molG)
                    BFG.add_node( frag, intensity=int(mol['estimate']) )

# how much of HTR should we find out? ANALYSIS OF THE HTR VS ETD
intensity_on_g = Counter()
for mols, error, status in MassTodonResults:
    if status=='optimal': #TODO what to do otherwise?
        for mol in mols:
            if mol['molType'][0] != 'p':
                intensity_on_g[mol['g']] += mol['estimate']
intensity_on_g
intensity_on_g[-1]/sum(intensity_on_g.values())
def protonate(Q,frag):
    a, b, c = {
        'p' : (1,0,1),
        'c' : (0,-1,0),
        'z' : (0,0,1)
    }[frag]
    for q in range(1,Q+a):
        for g in range(b,Q-q+c):
            yield (q,g)

cProt = set(protonate(Q,'c'))
zProt = set(protonate(Q,'z'))

cProt - zProt
zProt - cProt
zProt & cProt
###### Bloody 12 percent
