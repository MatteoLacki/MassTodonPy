%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules
from InSilico.spectrumGenerator import insilicoSpectrum, flatten, makeNoise
from operator import itemgetter
from PeakPicker.peakPicker import peakPicker
from PeakPicker.misc import deconvIter
from Solver.solver import solutions_iter

# getting random spectrum:
fasta = substanceP  = 'RPKPQQFFGLM'
Q = 3; modifications = {}; ionsNo = 1000000; P = .999
IS = insilicoSpectrum(fasta, Q, ionsNo, P)
sample,probs = IS.rvs(ionsNo)
MassSpectrum = [ (m, float(i))for m, i in flatten(sample)]
Noise = makeNoise( MassSpectrum, percentPeaks = .2 )

MassSpectrum.extend(Noise)
MassSpectrum.sort(key=itemgetter(0))

# peak picking like champions
chebCoverage = .99
jointProb = .999
precisionDigits = 3
precisionMass = .05 # In Daltons; the radius of mass

molecules_iter = genMolecules(fasta, 3, 'cz',modifications)
PP = peakPicker(    molecules_iter,
                    MassSpectrum,
                    chebCoverage,
                    jointProb,
                    precisionDigits )

# PP.add_M_n_E()
# PP.add_I()
# PP.add_eG()

G = PP.pickPeaks()
deconvProbs = deconvIter(G)
solutions   = solutions_iter(deconvProbs)

sols = list(solutions)

sols[0][0]

type(sols[0])
sols[0][1]
