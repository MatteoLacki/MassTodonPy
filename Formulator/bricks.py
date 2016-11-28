import igraph as ig
from linearCounter import LinearCounter as lCnt
from aminoAcid import AminoAcids

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
