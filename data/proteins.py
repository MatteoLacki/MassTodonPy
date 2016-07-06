import cPickle as pickle
from collections import Counter
from IsoSpecPy import IsoSpecPy
from math import exp
import numpy as np
import timeit


substanceP  = 'RPKPQQFFGLM'
ubiquitin   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

AA = pickle.load(open('data/AA.txt', 'rb'))
masses, probs = pickle.load(open('data/isotopes.txt', 'rb'))


def fasta2atoms(seq):   
	"""Get a counter of atoms from a protein fasta sequence."""
	ion = Counter()
	for i,s in enumerate(seq):
		ion += AA[s]    # += returns a new object, addition not in place
		if i==0:
			ion['H'] += 1
	ion['H'] += 1
	ion['O'] += 1
	return ion

# %time fasta2atoms('AACC')   
atoms = fasta2atoms(ubiquitin)
# print atoms

def getShiftedPrecursors(fasta, N):
	precursors = []
	for i in xrange(N):
		precursor = fasta2atoms(fasta)
		precursor['H'] += i
		precursors.append(precursor)
	return precursors


def getProbs(atomCnt):
  """Make isotopic spectrum for a given atom count."""
  elements = atomCnt.keys()    
  counts = [atomCnt[e] for e in elements]
  M = [ masses[e] for e in elements ]
  P = [ probs[e]  for e in elements ]
  res = IsoSpecPy.IsoSpec( counts, M, P, .95)
  envelope = [ (mass, exp(logProb)) for mass, logProb, _ in res.getConfs()]
  return envelope

E = getProbs(atoms)
print E

def randomSpectrum(E, sampleNo, massNoise=None, digit=2):
  """Simulate a spectrum: the mass of each ion is gaussian around the true mass.  
	
  Keyword arguments:
  E           -- the spectrum: a list of tuples (mass, probability)
  sampleNo    -- number of ions
  massNoise   -- the amount of noise around masses: spectrum resolution
  digit       -- the accuracy of mass representation: number of digits after coma
  """
  P = [p for m,p in E]
  M = [m for m,p in E]
  C = np.random.multinomial( sampleNo, P )
  NC = Counter()
  if massNoise:
      for m, c in zip(M,C): 
          for i in xrange(c):
              mass = round( np.random.normal(m, massNoise), digit)
              NC[mass] += 1
  else: 
      for m,c in zip(M,C):
          NC[round(m, digit)] += c
  return NC

sampleNo = 10000
massNoise = .3
randomSpectrum(E, sampleNo, massNoise, digit=2)


def randomlyCalibratedSpectrum(E, size, massNoise=None):
    """Generate a perfect spectrum subject only to mass shifts."""
    counts = np.random.multinomial( size, [prob for mass, prob in E] )
    if massNoise:
        masses = [mass+np.random.normal(0, massNoise) for mass, prob in E]
    else:
        masses = [mass for mass, prob in E]
    return zip(masses, counts)


def roundSpectrum(spec, digit):
    """Round spectrum masses to a given precision level."""
    M = Counter()
    for m, p in spec:
        M[ round(m,digit) ] += p
    masses = M.keys()
    masses.sort()
    return [(m, M[m]) for m in masses ]


def deconvolutionProblem(protein, numbersOfIons, massNoise=0.01, digit=2):
    """Prepare a task for MassTodon deconvolution."""
    N = len(numbersOfIons)
    envelopes = [getProbs(atomCnt) for atomCnt in getShiftedPrecursors(protein, N)]
    S = [ randomSpectrum(e, n, massNoise, digit) for e, n in zip(envelopes, numbersOfIons) ]
    empiria = Counter()
    for s in S:
      empiria += s
    M = list(empiria)
    M.sort()
    empiria = [ (m, empiria[m]) for m in M ]
    theory  = [ roundSpectrum(e,digit) for e in envelopes ]
    return (empiria, theory)


# numbersOfIons = [100000, 200000, 150000]
# empiria, theory = deconvolutionProblem(ubiquitin, numbersOfIons)


# import matplotlib.pyplot as plt
# import numpy as np

# mass, intensity = np.array(empiria).T
# plt.plot(mass, intensity)
# plt.show()

