from collections import defaultdict
from linearCounter.linearCounter import LinearCounter as lCnt

def standardize(modifications):
    '''Standardize modifications so that they meet the internal nomenclature scheme.'''
    backboneAtom2aaNomen = {'N':'L', 'Calpha':'C', 'C':'R'}
    R = defaultdict(lambda:defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        R[ tag[1]-1 ][ backboneAtom2aaNomen[tag[0]] ] = lCnt(atomCnt)
    return R

def countIsNegative(atomCnt):
    '''Check if any element of a dictionary is a negative number.'''
    return any( atomCnt[elem]<0 for elem in atomCnt )

def atomCnt2string(atomCnt):
    return "".join( el+str(cnt) for el, cnt in atomCnt.items() )
