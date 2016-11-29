%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules, makeFragments, pandizeSubstances
from Formulator.fasta2atomcnt import fasta2atomCnt
from InSilico.spectrumGenerator import insilicoSpectrum, genIsotopicEnvelope, flatten, makeNoise
from math import sqrt

fasta = substanceP  = 'RPKPQQFFGLM'
Q = 3; modifications = {}; ionsNo = 1000000; P = .999
molecules = list(genMolecules(fasta, 3, 'cz',modifications))
IS = insilicoSpectrum(fasta, Q, ionsNo, P)
sample, probs = IS.rvs(ionsNo)
MassSpectrum = [ (m, float(i))for m, i in flatten(sample)]
Noise = makeNoise(MassSpectrum, percentPeaks = .2)
MassSpectrum.extend(Noise)

from operator import itemgetter
MassSpectrum.sort(key=itemgetter(0))

import 	igraph as ig
import 	intervaltree as it
