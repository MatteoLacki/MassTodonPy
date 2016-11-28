import igraph as ig
import pandas as pd
from collections import Counter, defaultdict
from aminoAcid import AminoAcids
from linearCounter import LinearCounter as lCnt
from itertools import chain

def standardize(modifications):
    '''Standardize modifications so that they meet the internal nomenclature scheme.'''
    backboneAtom2aaNomen = {'N':'L', 'Calpha':'C', 'C':'R'}
    R = defaultdict(lambda:defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        R[ tag[1]-1 ][ backboneAtom2aaNomen[tag[0]] ] = lCnt(atomCnt)
    return R


def elementContent(G):
    '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
    atomNo = lCnt()
    for el in G.vs['elem']:
        atomNo[el] += 1
    return atomNo


def makeBricks():
    '''Prepare the base counters of atoms.

    These will be later collected into superatoms that finally make up the observed fragments.'''
    AAS = AminoAcids().get()
    AAS = [ (aa, AAS[aa]['graph'],AAS[aa]['NalphaIDX'],AAS[aa]['CcarboIDX']) for aa in AAS ]
    bricks = {}
    for aa, G, N, C in AAS:
        N = G.vs[N]['name']
        C = G.vs[C]['name']
        G.delete_edges(Roep_ne=None)
        G = G.decompose()
        brick = {}
        if len(G) == 2:
            assert aa=='P'
            brick['C'] = lCnt()
        while len(G) > 0:
            g = G.pop()
            isN = N in g.vs['name']
            isC = C in g.vs['name']
            if isN or isC:
                if isN:
                    tag = 'L'
                elif isC:
                    tag = 'R'
                else:
                    raise ValueError
            else:
                tag = 'C'
            brick[tag] = elementContent(g)
        bricks[aa] = brick
    return bricks


def countIsNegative(atomCnt):
    '''Check if any element of a dictionary is a negative number.'''
    return any( atomCnt[elem]<0 for elem in atomCnt )


def prolineBlockedFragments(fasta):
    '''Checks which c-z fragments cannot occur.'''
    blocked = set('c0')
    for i, f in enumerate(fasta):
        if f=='P':
            blocked.add( 'c' + str(i) )
            blocked.add( 'z' + str( len(fasta)-i ) )
    return blocked


def fasta2atomCnt(fasta, modifications = {}):
    '''Translates a fasta with modifications into a atom count.'''
    modifications = uniformify(modifications)
    bricks = makeBricks()
    def getBrick(aaPart):
        brick = bricks[aa][aaPart] + modifications[aaNo][aaPart]
        if countIsNegative(brick):
            print("Attention: your modification has an unexpected effect. Part of your molecule now has negative atom count. Bear that in mind while publishing your results.")
        return brick
    atomCnt = lCnt()
    for aaNo, aa in enumerate(fasta):
        atomCnt += getBrick('L') + getBrick('C') + getBrick('R')
    atomCnt += lCnt({'O':1,'H':2})
    return Counter(atomCnt)


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


def makeFragments(fasta, type, modifications):
    modifications = standardize(modifications)
    fragmentator = {
        'cz': make_cz_fragments
    }[type]
    return fragmentator(fasta, modifications)
