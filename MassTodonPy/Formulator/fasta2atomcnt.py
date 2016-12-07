# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from linearCounter import linearCounter as lCnt
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
