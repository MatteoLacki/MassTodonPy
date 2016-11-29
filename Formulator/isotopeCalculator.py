from IsoSpecPy import IsoSpecPy
from misc import atomCnt2string
from math import exp
from collections import Counter
try:
  import cPickle as pickle
except:
  import pickle

class IsotopeCalculations:
    def __init__(self, isoMasses=None, isoProbs=None):
        if isoMasses==None or isoProbs==None:
            isoMasses, isoProbs = pickle.load(open('data/isotopes.txt', 'rb'))
        self.isoMasses = isoMasses
        self.isoProbs  = isoProbs

        self.elementsMassMean = dict(
            (el, sum( pr*m for pr, m in zip(self.isoProbs[el], self.isoMasses[el]) ) )
            for el in self.isoMasses.keys() )

        self.elementsMassVar  = dict(
            (el, sum( pr*m**2 for pr, m in zip(self.isoProbs[el], self.isoMasses[el])) - self.elementsMassMean[el]**2 )
            for el in self.isoMasses.keys() )

    def getMonoisotopicMass(self, atomCnt):
        return sum( self.isoMasses[el][0]*elCnt for el, elCnt in atomCnt.items() )

    def getMassMean(self, atomCnt):
        return sum( self.elementsMassMean[el]*elCnt for el, elCnt in atomCnt.items() )

    def getMassVar(self, atomCnt):
        return sum( self.elementsMassVar[el]*elCnt for el, elCnt in atomCnt.items() )

    def getIsotopicEnvelope(self, atomCnt, jointProb, precisionDigits = 3):
        counts = []
        masses = []
        probs  = []
        for el, cnt in atomCnt.items():
            counts.append(cnt)
            masses.append(self.isoMasses[el])
            probs.append(self.isoProbs[el])
        I = IsoSpecPy.IsoSpec(tuple(counts), tuple(masses), tuple(probs), jointProb)
        result = Counter()
        for conf in I.getConfs():
            mz = conf[0]
            intensity = conf[1]
            result[ round(mz, precisionDigits) ] += exp(intensity)
        return list(result.items())
