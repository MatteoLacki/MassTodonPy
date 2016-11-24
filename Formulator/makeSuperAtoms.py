import igraph as ig
import pandas as pd
from collections import Counter, defaultdict
from aminoAcid import AminoAcids
from linearCounter import LinearCounter as lCnt

substanceP = 'RPKPQQFFGLM'
ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
testFasta = 'IKKLMNNMMM'

# fasta = substanceP
fasta = testFasta

backboneAtom2aaNomen = {'N':'L', 'Calpha':'C', 'C':'R'}
# modifications = {   ('N',2) :       lCnt({'H': -1, 'O': +2, 'N': +3}),
#                     ('Calpha',2) :  lCnt({'H': -1, 'O': +2, 'N': +3}),
#                     ('Calpha',5) :  lCnt({'H': -2, 'S': +2, 'N': +2}),
#                     ('C',6) :       lCnt({'H': -2, 'S': +2, 'N': +2}) }

# modifications = {   ('N',2) :       lCnt({'H': -1, 'O': +2, 'N': +3}),
                    # ('Calpha',2) :  lCnt({'H': -1, 'O': +2, 'N': +3}),
                    # ('Calpha',5) :  lCnt({'H': -2, 'S': +2, 'N': +2}),
                    # ('C',6) :       lCnt({'H': -2, 'S': +2, 'N': +2}) }
modifications = {}

modDiff = sum(modifications.values())

def uniformifyModifications(modifications):
    R = defaultdict(lambda:defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        R[ tag[1]-1 ][ backboneAtom2aaNomen[tag[0]] ] = atomCnt
    return R

modifications = uniformifyModifications(modifications)

def elementContent(G):
    '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
    atomNo = lCnt()
    for el in G.vs['elem']:
        atomNo[el] += 1
    return atomNo

def makeBricks():
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
    return any( atomCnt[elem]<0 for elem in atomCnt )

def make_cz_ions(fasta):
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

    precursor = sum(superAtoms)

    fragments = []
    N = len(superAtoms)
    cFrag = lCnt({'H':1}) # Adding one extra hydrogen to meet the definition of a c fragment.
    for i in range(N-1):
        cFrag += superAtoms[i]
        cFrag_tmp = lCnt(cFrag)
        cFrag_tmp['type'] = 'c'+str(i)
        fragments.append(cFrag_tmp)

    zFrag = lCnt()
    for i in range(1,N):
        zFrag += superAtoms[N-i]
        zFrag_tmp = lCnt(zFrag)
        zFrag_tmp['type'] = 'z'+str(i)
        fragments.append(zFrag_tmp)

    return precursor, fragments


def pandizeSubstances(precursor, fragments):
    precursor['type'] = 'precursor'
    result = fragments
    result.append(precursor)
    result = pd.DataFrame(result).fillna(0)
    result[list('CHNOS')] = result[list('CHNOS')].astype(int)
    idx = ['type']
    idx.extend(list('CHNOS'))
    result = result[idx]
    return result

pandizeSubstances(*make_cz_ions(fasta))
