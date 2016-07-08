import igraph as ig
import numpy  as np
import itertools as it
import pickle

from collections import Counter

aminoAcids  = pickle.load(open('aminoAcids.p', 'rb'))
ubiquitin   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
fasta       = ubiquitin

def elementContent(G):
    '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
    atomNo = Counter()  
    for el in G.vs['elem']:
        atomNo[el] += 1
    return atomNo

# def simplifyAminoAcid(aaTag, bonds2break = ('cz',)):
#     if isinstance( bonds2break, type(tuple()) ):
#         pass
#     elif isinstance( bonds2break, type(dict()) ):
#         pass
#     else:
#         print('Unrecognizable format of bonds to break.')
#         pass

def establishFragmentType(fragment, aaType): 
    '''Establish what is the fragment type: (L)eft, (C)enter, or (R)ight, or (LC) or (R) for proline.'''
    if aaType == 'P':
        if 'Ccarbo' in fragment.vs['name']:
            return 'R'
        else:
            return 'LC'
    else:
        if 'Calpha' in fragment.vs['name']: 
            return 'C'
        elif 'Nalpha' in fragment.vs['name']: 
            return 'L'
        else:
            return 'R'


def getSuperAtoms(fasta, fragmentTypes):
    '''Makes super atoms - basic bricks of differences that should be observed in the experimental spectrum.'''
    fragments = {}
    for f in set(fasta):
        G = aminoAcids[f]['graph'].copy()
        G.delete_edges(Roep_ne=None)        
        G = G.decompose()
        G = dict( (establishFragmentType(g,f), elementContent(g)) for g in G ) 
        if set(G) == set(['L','C','R']):
            if 'ax' in fragmentTypes: 
                if 'cz' in fragmentTypes:
                    G = [ ('LL', G['L']), ('CC',G['C']), ('RR',G['R']) ]
                else:
                    G = [ ('LC', G['L']+G['C']), ('RR', G['R']) ]
            else:
                if 'cz' in fragmentTypes:
                    G = [ ('LL', G['L']), ('CR', G['C']+G['R']) ]
                else:
                    G = [ ('LR', G['L']+G['C']+G['R']) ]    
        elif set(G) == set(['LC','R']): 
            if 'ax' in fragmentTypes:
                G = [ ('LC', G['LC']), ('RR', G['R']) ]            
            else:
                G = [ ('LR', G['LC']+G['R']) ]
        else:
            print('Fuck a doodle doo!')
        fragments[f] = G
    #
    superAtoms = [] # It must be a list to be mutable.
    for fNo, f in enumerate(fasta):
        if fNo==0 or 'by' in fragmentTypes:
            for aaType, atomCnt in fragments[f]:    
                superAtoms.append([ aaType, fNo, fNo, atomCnt.copy() ])
        else:
            for fragNo, fragment in enumerate(fragments[f]):
                aaType, atomCnt = fragment
                if fragNo==0:
                    # [ aaType, fNo, fNo, Counter() ]
                    superAtoms[-1][-1].update(atomCnt)              
                    # Update left right fragment tags.      
                    superAtoms[-1][0] = superAtoms[-1][0][0] + aaType[-1]
                    # Upadate second fragment tag counter.
                    superAtoms[-1][2] = fNo
                else:
                    superAtoms.append([ aaType, fNo, fNo, atomCnt.copy() ])
    #
    # Adding water to the molecule: acids' backbones were water free
    superAtoms[0][-1]['H']  += 1
    superAtoms[-1][-1]['H'] += 1
    superAtoms[-1][-1]['O'] += 1
    #
    assert any( sA[1]<=sA[2] for sA in superAtoms )
    return superAtoms

def makeFragments(fasta, fragmentTypes=['cz'], innerFragments = False):
    '''Makes tagged chemical formulas of fragments under given fragmentation scheme.'''    
    fragmentTypes   = set(fragmentTypes)
    superAtoms      = getSuperAtoms(fasta, fragmentTypes)
    fragments = []
    if innerFragments:
        for i in range(len(superAtoms)):
            fragments.append(superAtoms[i])
            for j in range(i+1,len(superAtoms)):
                SA1     = fragments[-1]
                SA2     = superAtoms[j]
                aaType  = SA1[0][0]+SA2[0][1]
                lFragNo = SA1[1]
                rFragNo = SA2[2]
                atomCnt = SA1[3]+SA2[3]
                fragments.append([aaType, lFragNo, rFragNo, atomCnt])
    else:
        prevCnt = Counter()
        # abc fragments--------------------------------------------------
        for aaType, lFragNo, rFragNo, atomCnt in superAtoms:  
            aaType  = 'L'+aaType[-1]
            lFragNo = 0
            prevCnt = prevCnt + atomCnt
            fragments.append([aaType, lFragNo, rFragNo, prevCnt])
        #----------------------------------------------------------------
        precursor1  = fragments[-1]
        prevCnt     = Counter()
        superAtoms.reverse()
        L = len(fasta)
        # xyz fragments--------------------------------------------------
        for aaType, lFragNo, rFragNo, atomCnt in superAtoms:  
            aaType  = aaType[0]+'R'
            rFragNo = L-1
            prevCnt = prevCnt + atomCnt
            fragments.append([aaType, lFragNo, rFragNo, prevCnt])
        #----------------------------------------------------------------
        precursor2  = fragments[-1]
        assert precursor1==precursor2, print('Precursors done from left and right differ.')
        # Removing extra precursor---------------------------------------
        del fragments[-1]
    return fragments

print(
len(makeFragments( fasta,['cz'] )),
len(makeFragments( fasta,['cz','ax'] )),
len(makeFragments( fasta,['cz','ax','by'],innerFragments=True)),
len(makeFragments( fasta,['cz',],innerFragments=True)),
len(makeFragments( fasta,['cz','ax',],innerFragments=True))
)

# Check this bloody function.
def roepstorffy(fragment,fasta):
    '''Convert my naming convention into something similar to Roepstorff fragmentation conventation.'''
    lcr2abc = {'L':'c','C':'a','R':'b'}
    lcr2xyz = {'L':'y','C':'z','R':'x'}
    L = len(fasta)
    AAtype, AAleft, AAright, atomCnt = fragment
    lType, rType = AAtype
    if lType=='L' and AAleft==0:
        nameL = ''
    else:
        lType = lcr2xyz[lType]
        nameL = lType + str(L-AAleft+(0 if lType=='x' else 1))    
    if rType=='R' and AAleft==L-1:
        nameR = ''
    else:
        rType = lcr2abc[rType]
        nameR = rType + str(AAright-(1 if rType=='c' else 0))
    return nameL+nameR 
# roepstorffy(superAtoms[10], fasta)

