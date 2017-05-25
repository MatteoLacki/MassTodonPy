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
from MatchMaker         import czMatchMakerBasic as analyzer_basic, czMatchMakerIntermediate as analyzer_inter, czMatchMakerUpperIntermediate as analyzer_up_inter
from collections        import Counter

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
        self.modifications = modifications


    def readSpectrum(   self,
                        spectrum=None,
                        path=None,
                        cutOff=100.,
                        digits=2,
                        topPercent=1.0 ):
        '''Read in a mass spectrum.

        Read either an individual text file or merge runs from an mzXml files. In case of the mzXml file
        '''
        self.spectrum, self.total_spectrum_intensity = readSpectrum(path, spectrum, cutOff, digits, topPercent)


    def prepare_problems(self, M_minProb=.75):
        '''Prepare a generator of deconvolution problems.'''
        self.problems = self.peakPicker.get_problems(self.spectrum, M_minProb)

        #TODO: add multiprocessing: turn of BLAS asynchronic calculations
        #TODO: make a more sensible use of sequentiality
    def run(self, solver='sequential', method='MSE', max_times_solve=5, **args):
        '''Perform the deconvolution of problems.'''
        self.res = solve(problemsGenerator = self.problems,
                    args   = args,
                    solver = solver,
                    method = method,
                    max_times_solve = max_times_solve )


    def summarize_results(self):
        summary = Counter()
        for r in self.res:
            summary['L1_error'] += r['L1_error']
            summary['L2_error'] += r['L2_error']
            summary['underestimates']+= r['underestimates']
            summary['overestimates'] += r['overestimates']
        summary['relative_L1_error'] = summary['L1_error']/self.total_spectrum_intensity
        summary['relative_L2_error'] = summary['L2_error']/self.total_spectrum_intensity
        summary['%% percent of explained spectrum'] = 1.0 - summary['underestimates']/self.total_spectrum_intensity
        return summary


    def flatten_results(self, minimal_estimated_intensity=100.0):
        '''Return one list of results, one list of difficult cases, and the error.'''
        optimal     = []
        nonoptimal  = []
        totalError  = 0.0
        for mols, error, status in self.res:
            if status=='optimal':
                totalError += error
                for mol in mols:
                    if mol['estimate'] > minimal_estimated_intensity:
                        mol_res = {}
                        for key in ['estimate', 'molType', 'q', 'g', 'formula']:
                            mol_res[key] = mol[key]
                        optimal.append(mol_res)
            else:
                nonoptimal.append(mols)
        return optimal, nonoptimal, totalError


    #TODO: push Ciach to write a general structure.
    def get_subsequence(self, name):
        '''For purpose of ETDetective, to bridge gaps in definitions.'''
        if self.modifications == {'C11': {'H': 1, 'N': 1, 'O': -1}}:
            suffix = '*'
        else:
            suffix = ''
        if name[0]=='p':
            return '*'+self.fasta+suffix
        if name[0]=='z':
            return self.fasta[ len(self.fasta)-int(name[1:]): ] + suffix
        else:
            return '*' + self.fasta[ 0:int(name[1:]) ]


    def i_flatten_results_for_ETDetective(self):
        '''Generate pairs substance-estimate.'''
        for subproblem in self.res:
            for r in subproblem['alphas']:
                seq = self.get_subsequence( r['molType'] )
                yield (seq, r['q'], r['g'], r['molType']), r['estimate']


    def gen_ETDetective_inputs(self):
        '''Get inputs for ETDetective.'''
        data = {}
        for r in self.i_flatten_results_for_ETDetective():
            (seq, q, g, molType), estimate = r
            if q == self.Q and g == 0 and molType=='precursor':
                precursor = seq, q, g, molType
            else:
                data[(seq, q, g, molType)] = estimate
        return precursor, data


    def analyze_reactions(self, analyzer='basic', accept_nonOptimalDeconv = False, min_acceptEstimIntensity = 100., verbose=False, **advanced_args):
        '''Estimate reaction constants and quantities of fragments.'''
        chosen_analyzer = {
            'basic': analyzer_basic,
            'inter': analyzer_inter,
            'up_inter': analyzer_up_inter
        }[analyzer](self.res, self.Q, self.fasta,
                    accept_nonOptimalDeconv,
                    min_acceptEstimIntensity, verbose )
        return chosen_analyzer.pair()
