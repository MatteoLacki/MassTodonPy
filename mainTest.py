%load_ext autoreload
%autoreload
from Formulator.formulator import makeFragments, pandizeSubstances, fasta2atomCnt
from Formulator.protonations import protonate
from itertools import chain
from IsoSpecPy import IsoSpecPy
from math import exp, sqrt
from Formulator.isotopeCalculator import IsotopeCalculations

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

precursor, cFrags, zFrags = makeFragments(fasta, 'cz', modifications)
# pandizeSubstances(precursor, cFrags, zFrags)



# def getMonoisotopicMass(mol):

mol = list(precursor())[0]




IC = IsotopeCalculations()
IC.getMonoisotopicMass(mol['atomCnt'])
IC.getMassMean(mol['atomCnt'])
IC.getMassVar(mol['atomCnt'])

for mol in chain(precursor(),cFrags(),zFrags()):
    for q,g in protonate(Q,mol['type']):
        print q,g,mol


i = IsoSpecPy.IsoSpec.IsoFromFormula("H2O1", 0.9)

defProtonizeMolecules

def C():
    for c in cFrags():
        for q,c in protonate(Q,'c'):
            print c
            c['q'] = q
            c['g'] = g
            yield c


list(protonate(Q,'precursor'))



print i.getConfs()
