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
# .7O......O............................ZZZZZZZ...................................
# .7......I$...............................Z......................................
# .7.Z....O7....,ZZZ.....OZZ......ZZO .....Z.....OZZ......OZO.....,ZZO....ZOOZ+...
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

    def randomSpectrum(
            self,
            ionsNo,
            aaPerOneCharge  = 5,
            jointProb       = .999,
            scale           = .01,
            percentPeaks    = .2     ):

        masses, intensities = self.isoCalc.randomFragmentationExperiment( self.fasta, self.Q, ionsNo, self.formulator, aaPerOneCharge, jointProb, scale )

        noise_masses, noise_intensities = self.isoCalc.addNoise( masses, intensities, percentPeaks )

        return masses, intensities, noise_masses, noise_intensities

        #TODO: add multiprocessing
    def run(self,
            spectrum,
            M_minProb   = .75,
            solver      = 'sequential',
            method      = 'MSE',
            **args
        ):
        problems = self.peakPicker.get_problems(spectrum, M_minProb)
        return solve(problems,solver,method,**args)
