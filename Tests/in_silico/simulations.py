%load_ext autoreload
%load_ext line_profiler
%autoreload
from MassTodonPy import MassTodon
from MassTodonPy.Formulator import makeFormulas
from collections import Counter

fasta = 'RPKPQQFFGLM'
Q = 3; jP= .999; mzPrec = .05; precDigits = 2; L2 = 0.00001; M_minProb = .7

M = MassTodon(  fasta = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

Formulator = makeFormulas(fasta, Q)
mols = list(Formulator.makeMolecules())
atomCnts = dict( (molType, atomCnt) for molType, atomCnt, _, q, g in mols)
mols = set( (molType, atomCnt, q, g) for molType, atomCnt, _, q, g in mols)

observedMols = Counter()

precursor = ('precursor', atomCnts['precursor'], )
atomCnts['precursor']
