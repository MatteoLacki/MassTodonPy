#
# -*- coding: utf-8 -*-
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
# ................................................................................
# ............~MMM:...............................................................
# ..........:M....NM..............................................................
# ..........8M,....M,.............................................................
# ............M~...8O.............................................................
# ..........+.MM,..,M..............................?MMMMMN=.......................
# ..........:..:M:..M..........................DMMM~..  ....MMMM$.................
# ..........+...DN..M...............+MMMMMMMMMMM............... $MMM=.............
# ..........?...OM..M=............ZMD...... .......................ZMM?...........
# ..........:...?M..8O...........M,...................................MM,.........
# ..........+...7M..~M..........M......................................OMM........
# ..........+...IM...M.........MI....................................... MM.......
# ..........:...ZM...M.........M...M?.8D..................................MM......
# ..........:...DM...M:.......M:..M....,O.................................?M:.....
# ..............NM...~M......OM...M....8:..................................MM.....
# ...........$..ZM....MO.....M.....DMMM.....................................M,....
# ...........Z..,M.....MM..8M:..............................................MO....
# ..........:Z...MN......DMO................................................MM....
# ..........:Z...~M~........................................................7M....
# ..........?Z....$M7........................................................M....
# ..........+Z.....=MM...........MZ..........................................,M...
# ..........~Z.......MMMM=....IMMMM..~M....................................:M.M...
# ..........=Z......:...IMMMMMZ..MI..NMM....................................NM=...
# ..........+Z......MMZ.........OM...M,.M....................................M+...
# ..........:Z.....Z.M.,MMMIIDMMM....M..M.N~..................................M...
# ..........:Z.....Z..MM...........~MD..M.MM................................M,M7..
# ..........:Z.....Z....7MMMO:~8MMM7...ND.MN................................MM:+..
# ..........:Z.....Z..................?I..MD................................M.M...
# ..........=$.....?.....................MZN................................M+....
# ..........ZZ.....7....................MD.M........MMM..N,..M...M7M........MO....
# ..........Z$~....?.......................M........M..M.~D.~N..:,.M........MO....
# ..........ZZ,....$.... .......=..........M........M:.M..D..M..+=.M,.......MO....
# ........~.ZZ?....$....=.......=..........MM......?M,.,...........MM......8M=....
# .$I.,.IZ,I.+,~I.,....~~..7Z.Z.:...=.,$Z...~8MMMMMMI..............~MMMMMMMM8.....
# ................................................................................
# ................................................................................
# ................................................................................
# .7O......O............................ZZZZZZZ..............Z....................
# .7......I$...............................Z.................Z....................
# .7.Z....O7....,ZZZ.....OZZ......ZZO .....Z.....OZZ......OZOZ....,ZZO....ZOOZ+...
# .7..$..Z.7...7....Z...Z....=.......O.....Z....O....?...=...?...I....Z...Z....O..
# .7..O..~.7........Z.. Z........7.........Z....O....O...O...O...Z....Z...O....Z..
# .7...ZZ..7.....8Z.Z.....O.......$Z.......Z....O....O...O...O...Z....Z...O....Z..
# .7...~...7...Z....Z.......7........O.....Z....O....O...O...O...Z....Z...O....Z..
# .7.......7...O....O...O... I.......Z.....Z....Z... $...O...O...$....O...O....Z..
# .7.......7....OOOOO....ZOZI....:ZOZ......Z.....ZOZ7.....OZZ.....7ZOZ....Z....Z..
# ................................................................................

from Formulator         import makeFormulas
from IsotopeCalculator  import isotopeCalculator
from PeakPicker         import PeakPicker
from Solver             import solve
from Parsers            import readSpectrum

class MassTodon():
    def __init__(   self,
                    fasta,
                    precursorCharge,
                    fragType        = 'cz',
                    precDigits      = 2,
                    mzPrec          = .05,
                    jointProbability= 0.999,
                    isoMasses       = None,
                    isoProbs        = None,
                    modifications   = {} ):
        self.fasta  = fasta
        self.Q      = precursorCharge

        self.Forms  = makeFormulas(
            fasta   = fasta,
            Q       = self.Q,
            fragType= fragType,
            modifications = modifications )

        self.IsoCalc = isotopeCalculator(
            jP          = jointProbability,
            precDigits  = precDigits,
            isoMasses   = isoMasses,
            isoProbs    = isoProbs )

        self.peakPicker = PeakPicker(
            Forms   = self.Forms,
            IsoCalc = self.IsoCalc,
            mzPrec  = mzPrec )


    def readSpectrum(   self,
                        spectrum=None,
                        path=None,
                        cutOff=100.,
                        digits=2,
                        topPercent=1.0 ):
        '''Read in a mass spectrum.

        Read either an individual text file or merge runs from an mzXml files. In case of the mzXml file
        '''
        if path:
            self.spectrum = readSpectrum( path,
            cutOff, digits, topPercent)
        else:
            if spectrum:
                self.spectrum = spectrum # a tupple of two numpy arrays
            else:
                print 'Wrong path or you did not provide a spectrum.'
                raise KeyError

    def prepare_problems(self, M_minProb=.75):
        '''Prepare a generator of deconvolution problems.'''
        self.problems = self.peakPicker.get_problems(self.spectrum, M_minProb)

        #TODO: add multiprocessing
    def run(self, solver='sequential', method='MSE', **args):
        '''Perform the deconvolution of problems.'''
        res = solve(
            problemsGenerator = self.problems,
            solver = solver,
            method = method,
            args = args)
        return res
