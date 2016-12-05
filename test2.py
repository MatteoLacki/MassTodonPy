%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules
from InSilico.spectrumGenerator import insilicoSpectrum, flatten, makeNoise
from operator import itemgetter
from PeakPicker.peakPicker import peakPicker
from PeakPicker.misc import deconvIter
from Solver.solver import solutions_iter
try:
   import cPickle as pickle
except:
   import pickle

# getting random spectrum:
fasta = substanceP = 'RPKPQQFFGLM'
fasta = ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9

modifications = {}
ionsNo  = 1000000
P       = .999

IS = insilicoSpectrum(fasta, Q, ionsNo, P)

f_charges       = [ float(mol[1]**2) for mol in IS.molecules ]
total_f_charges = sum(f_charges)
probs           = [ q/total_f_charges for q in f_charges ]

from numpy.random import multinomial
moleculeCounts = multinomial( ionsNo, probs )

from IsoSpecPy import IsoSpecPy

def genIsotopicEnvelope(isotopologuesNo, atomCnt, jointProb=.999):
    atomCnt_str = atomCnt2string(atomCnt)
    envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
    return

def atomCnt2string(atomCnt):
    return "".join( el+str(cnt) for el, cnt in atomCnt.items() )

%time sample = [ genIsotopicEnvelope(isoCnt, mol[3], IS.P) for isoCnt, mol in zip(moleculeCounts, IS.molecules)]


from collections import Counter

x = Counter({2.4:3, 2.5:6})
y = Counter({2.6:3, 2.5:6})




import multiprocessing as mu
def getEnv(args):
    return genIsotopicEnvelope( args[0], args[1], .999)

p = mu.Pool(3)
p.map(getEnv, zip(moleculeCounts, IS.molecules))



sample,probs = IS.rvs(ionsNo)
MassSpectrum = [ (m, float(i))for m, i in flatten(sample)]
Noise = makeNoise( MassSpectrum, percentPeaks = .2 )

with open('data/spectrum.spec','w') as f:
    pickle.dump( (MassSpectrum, Noise), f)

MassSpectrum.extend(Noise)
MassSpectrum.sort(key=itemgetter(0))



# peak picking like champions
chebCoverage    = .99
jointProb       = .999
precisionDigits =  3
precisionMass   = .05 # In Daltons; the radius of mass

molecules_iter = genMolecules(fasta, Q, 'cz',modifications)
X = list(molecules_iter)
X



PP = peakPicker(    molecules_iter,
                    MassSpectrum,
                    chebCoverage,
                    jointProb,
                    precisionDigits )

G = PP.pickPeaks()
deconvProbs = deconvIter(G)
solutions   = solutions_iter(deconvProbs)
alphaType   = frozenset('MR')
alphas = []
explainedRatios = []

successful      = []
not_successful  = []
for sol in solutions:
    res, g, coord2edges, linprogInput = sol
    if res.success:
        successful.append(sol)
        explainedRatios.append( -res.fun/( res.slack.sum()-res.fun ) )
        for k, (n1, n2) in enumerate(coord2edges):
            if frozenset((n1[0], n2[0])) == alphaType:
                if n1[0]=='M':
                    M = g.node[n1]
                else:
                    M = G.node[n2]
                M['est_intensity'] = res.x[k]
                alphas.append(M)
    else:
        not_successful.append(sol)

len(alphas)
explainedRatios
