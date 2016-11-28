%load_ext autoreload
%autoreload
# from Formulator.formulator import makeFragments, pandizeSubstances, fasta2atomCnt
# from Formulator.protonations import protonate
# from itertools import chain
# from IsoSpecPy import IsoSpecPy
# from math import exp, sqrt
# from Formulator.isotopeCalculator import IsotopeCalculations

from Formulator.formulator import genMolecules

ubiquitin   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
substanceP  = 'RPKPQQFFGLM'
fastas      = [substanceP, ubiquitin, ubiquitin+ubiquitin+ubiquitin]
fasta = substanceP
Q = 3
# modifications = {   ('N',2) :       {'H': -1, 'O': +2, 'N': +3},
#                     ('Calpha',2) :  {'H': -1, 'O': +2, 'N': +3},
#                     ('Calpha',5) :  {'H': -2, 'S': +2, 'N': +2},
#                     ('C',6) :       {'H': -2, 'S': +2, 'N': +2} }
modifications = {}
# pandizeSubstances(precursor, cFrags, zFrags)

molecules = list( genMolecules(fasta, 3, 'cz', modifications) )

from numpy.random import multinomial


ionsNo = 100000
f_charges = [ float(mol[1]**2) for mol in molecules ]
total_f_charges = sum(f_charges)
probs  = [ q/total_f_charges for q in f_charges ]

moleculeNo = multinomial( ionsNo, probs )

def atomCnt2string(atomCnt):
    return "".join( el+str(cnt) for el, cnt in atomCnt.items() )

atomCnt = molecules[0][3]
jointProb =.999
N = moleculeNo[0]

def simulateIsotopicEnvelope(isotopologuesN, atomCnt, jointProb=.999):
    atomCnt_str = atomCnt2string(atomCnt)
    envelope = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
    masses = []
    probs = []
    for x in envelope.getConfs():
        masses.append(x[0])
        probs.append(exp(x[1]))
    counts = multinomial( isotopologuesN, probs )
    return list(zip(masses, counts))

simulateIsotopicEnvelope(N, atomCnt)
