%load_ext autoreload
%load_ext line_profiler
%autoreload
from IsotopeCalculator  import isotopeCalculator, formulaParser
from IsoSpecPy          import IsoSpecPy
from math               import exp, floor, fsum

Q = 8; modifications = {}; jointProb = .999; mzPrec = .05; massPrecDigits = 2
isoMasses, isoProbs = None, None
isoCalc = isotopeCalculator(0, isoMasses, isoProbs)

atomCnt_str = 'C100H200'
isoCalc.isoEnvelope(atomCnt_str, jointProb, 3, 1)

masses, probs = isoCalc.getEnvelope(atomCnt_str, jointProb)


g = 1
q = 3
masses2 = np.around( (masses + g + q)/q, decimals=0 )


def aggregate2( keys, values, digits=2 ):
    '''Aggregate values with the same keys.'''
    uniqueKeys, indices = np.unique( keys, return_inverse=True)
    return uniqueKeys, np.bincount( indices, weights=values )


aggregate2(masses2, probs, digits=0)





isotope_masses = []
counts = []
isotope_probs  = []

formParser  = formulaParser()
atomCnt     = formParser.parse(atomCnt_str)

for el, cnt in atomCnt.items():
    counts.append(cnt)
    isotope_masses.append(isoCalc.isoMasses[el])
    isotope_probs.append(isoCalc.isoProbs[el])

envelope = IsoSpecPy.IsoSpec( counts, isotope_masses, isotope_probs, jointProb )
masses, logprobs, _ = envelope.getConfsRaw()

def cdata2numpyarray(x):
    '''Turn c-data into a numpy array.'''
    res = np.empty(len(x))
    for i in xrange(len(x)):
        res[i] = x[i]
    return res

masses  = cdata2numpyarray(masses)
probs   = np.exp(cdata2numpyarray(logprobs))

def agg_spec_proper(masses, probs, digits=2):
    '''Aggregate values with the same keys.'''
    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array(lists.keys())
    newProbs  = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs,op_flags=['readwrite']), newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs

agg_spec_proper(masses, probs, 0)



# memoization
self.isotopicEnvelopes[ atomCnt_str ] = ( masses, probs )
return masses.copy(), probs.copy()
