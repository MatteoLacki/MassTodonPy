import igraph as ig
from collections import Counter
from aminoAcid import AminoAcids
try:
  import cPickle as pickle
except:
  import pickle


fasta = ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
modifications = {   ('N', 10) :     Counter({'H': -1, 'O': +2, 'N': +3}),
                    ('Calpha',52) : Counter({'H': -2, 'S': +2, 'N': +2}),
                    ('C',69) :      Counter({'H': -2, 'S': +2, 'N': +2}) }

def elementContent(G):
    '''Extracts numbes of atoms of elements that make up the graph of a molecule.'''
    atomNo = Counter()
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
            brick['C'] = Counter()
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

bricks = makeBricks()




# def make_cz_ions(fasta):
superAtoms = []
sA = Counter()
for aaNo, aa in enumerate(fasta):
    sA += bricks[aa]['L']+bricks[aa]['C']
    superAtoms.append(sA)
    sA = Counter(bricks[aa]['R'])
sA += Counter({'O':1,'H':1})
superAtoms.append(sA)
superAtoms[0] += Counter({'H':1})




modifications

G.delete_edges
G.es['Roep']
G.delete_edges(Roep_ne=None)

G.decompose()
