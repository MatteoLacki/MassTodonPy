import pandas as pd
from linearCounter import LinearCounter as lCnt
from itertools import chain
from isotopeCalculator import IsotopeCalculations
from protonations import protonate
from bricks import makeBricks
from misc import standardize, countIsNegative

def prolineBlockedFragments(fasta):
    '''Checks which c-z fragments cannot occur.'''
    blocked = set('c0')
    for i, f in enumerate(fasta):
        if f=='P':
            blocked.add( 'c' + str(i) )
            blocked.add( 'z' + str( len(fasta)-i ) )
    return blocked

def make_cz_fragments(fasta, modifications):
    '''Prepares the precursor and the c and z fragments atom counts.'''
    bricks = makeBricks()
    def getBrick(aaPart):
        brick = bricks[aa][aaPart] + modifications[aaNo][aaPart]
        if countIsNegative(brick):
            print("Attention: your modification has an unexpected effect. Part of your molecule now has negative atom count. Bear that in mind while publishing your results.")
        return brick

    superAtoms = []
    sA = lCnt()
    for aaNo, aa in enumerate(fasta):
        sA += getBrick('L')
        superAtoms.append( sA )
        sA = getBrick('C') + getBrick('R')
    sA += lCnt({'O':1,'H':1})
    superAtoms.append(sA)
    superAtoms[0] += lCnt({'H':1})
    N = len(superAtoms)

    def getPrecursor():
        precursor = sum(superAtoms)
        yield {'moleculeType': 'precursor', 'atomCnt': dict(precursor), 'sideChainsNo' : len(fasta), 'type':'p' }

    blockedFragments = prolineBlockedFragments(fasta)

    def getCfrags():
        cFrag = lCnt({'H':1}) # Adding one extra hydrogen to meet the definition of a c fragment.
        for i in range(N-1):
            cFrag += superAtoms[i]
            cFrag_tmp = lCnt(cFrag)
            fragType = 'c'+str(i)
            if not fragType in blockedFragments and not i == 0:
                yield {'moleculeType': fragType, 'atomCnt': dict(cFrag_tmp), 'sideChainsNo' : i, 'type':'c' }

    def getZfrags():
        zFrag = lCnt()
        for i in range(1,N):
            zFrag += superAtoms[N-i]
            zFrag_tmp = lCnt(zFrag)
            fragType = 'z'+str(i)
            if not fragType in blockedFragments:
                yield {'moleculeType': fragType, 'atomCnt': dict(zFrag_tmp), 'sideChainsNo' : i, 'type':'z' }

    return getPrecursor, getCfrags, getZfrags


def pandizeSubstances(precursor, cFrags, zFrags):
    '''Turns results into a pandas data frame.'''
    def combineMolecules():
        for x in chain(precursor(), cFrags(), zFrags()):
            x['atomCnt']['moleculeType'] = x['moleculeType']
            yield x['atomCnt']
    result = pd.DataFrame(combineMolecules()).fillna(0)
    result[list('CHNOS')] = result[list('CHNOS')].astype(int)
    idx = ['moleculeType']
    idx.extend(list('CHNOS'))
    result = result[idx]
    return result


def makeFragments(fasta, type='cz', modifications={}):
    '''Generate all possible fragments given a Roepstorf Scheme.
    '''
    modifications = standardize(modifications)
    fragmentator = {
        'cz': make_cz_fragments
    }[type]
    return fragmentator(fasta, modifications)


def genMolecules(fasta, Q, fragmentationScheme='cz', modifications={}, aaPerOneCharge= 5):
    '''Generate protonated molecules following a given fragmentation scheme.
    '''
    IC = IsotopeCalculations()
    precursor, cFrags, zFrags = makeFragments(fasta, fragmentationScheme, modifications)
    for mol in chain(precursor(),cFrags(),zFrags()):
        for q,g in protonate( Q, mol['type'] ):
            if q * aaPerOneCharge < mol['sideChainsNo']:
                atomCnt = dict(mol['atomCnt'])
                atomCnt['H'] += q + g
                monoisotopicMass= IC.getMonoisotopicMass(atomCnt)/float(q)
                massMean = IC.getMassMean(atomCnt)/float(q)
                massVar  = IC.getMassVar(atomCnt)/float(q**2)
                yield ( mol['moleculeType'], q, g, atomCnt, monoisotopicMass, massMean, massVar )
