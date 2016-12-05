from linearCounter.linearCounter import LinearCounter as lCnt
from aminoAcid import AminoAcids
from misc import standardize, countIsNegative
from bricks import makeBricks

def fasta2atomCnt(fasta, modifications = {}):
    '''Translates a fasta with modifications into a atom count.'''
    modifications = standardize(modifications)
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
    return dict(atomCnt)
