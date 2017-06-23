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

from Formulator         import make_formulas
from IsotopeCalculator  import IsotopeCalculator
from PeakPicker         import PeakPicker
from Solver             import solve
from Parsers            import read_n_preprocess_spectrum
from MatchMaker         import czMatchMakerBasic as analyzer_basic, czMatchMakerIntermediate as analyzer_intermediate, czMatchMakerAdvanced as analyzer_advanced
from Visualization      import ResultsPlotter
from Summarator         import summarize_results
from itertools          import izip
from math               import ceil, log10
from intervaltree       import Interval as interval, IntervalTree

class MassTodon():
    def __init__(   self,
                    fasta,
                    precursor_charge,
                    mz_prec,
                    modifications   = {},
                    frag_type       = 'cz',
                    joint_probability_of_envelope = 0.999,
                    iso_masses      = None,
                    iso_probs       = None,
                    verbose         = False
        ):
        '''Make MassTodon somewhat less extinct by creating its instance.'''

        self.mz_prec = mz_prec
            # precision one digit lower than the input precision of spectra, eg.
            # mz_prec = .05     -->     prec_digits = 3
            # mz_prec = .005    -->     prec_digits = 4
        # self.prec_digits = int(ceil(-log10(mz_prec)))+1
        self.prec_digits = int(ceil(-log10(mz_prec)))
        self.fasta  = fasta
        self.Q      = precursor_charge
        self.verbose= verbose

        self.Forms  = make_formulas(
            fasta   = fasta,
            Q       = self.Q,
            frag_type     = frag_type,
            modifications = modifications )

        self.IsoCalc = IsotopeCalculator(
            jP          = joint_probability_of_envelope,
            prec_digits = self.prec_digits,
            iso_masses  = iso_masses,
            iso_probs   = iso_probs,
            verbose     = verbose   )

        self.peakPicker = PeakPicker(
            _Forms   = self.Forms,
            _IsoCalc = self.IsoCalc,
            mz_prec = mz_prec,
            verbose = verbose   )

        self.ResPlotter = ResultsPlotter(mz_prec)
        self.modifications = modifications
        self.spectra = {}


    def read_n_preprocess_spectrum(self,
            path    = None,
            spectrum= None,
            cut_off = None,
            opt_P   = None
        ):
        '''Read in a mass spectrum and round the masses to prec_digits digits after 0.

        Read either an individual text file or merge runs from an mzXml files. In case of the mzXml file
        '''
        self.spectra= read_n_preprocess_spectrum(
            path    = path,
            spectrum= spectrum,
            prec_digits = self.prec_digits,
            cut_off = cut_off,
            opt_P   = opt_P      )

        if self.verbose:
            print
            print 'original total intensity',   self.spectra['original total intensity']
            print 'total intensity after trim', self.spectra['total intensity after trim']
            print 'trimmed intensity', self.spectra['trimmed intensity']
            print

    # TODO is the thing below necessary?
    def spectrum_iter(self, spectrum_type):
        assert spectrum_type in ['original', 'trimmed'], "No such kind of spectrum: %s." % spectrum_type
        for mz, intensity in izip(*self.spectra[spectrum_type]):
            yield {'mz':mz, 'intensity':intensity}


    def run(self,
            solver              = 'sequential',
            multiprocesses_No   = None,
            method              ='MSE',
            max_times_solve     = 5,
            min_prob_per_molecule = .75,
            bootstrap           = False,
            forPlot             = False,
            **args ):
        '''Perform the deconvolution of problems.'''


        picker_res = self.peakPicker.get_problems(
            massSpectrum            = self.spectra['trimmed'],
            min_prob_per_molecule   = min_prob_per_molecule,
            bootstrap               = bootstrap )
        if bootstrap:
            self.problems, self.small_graphs_no_G = picker_res
        else:
            self.problems = picker_res

        self.res = solve(   problems = self.problems,
                            args   = args,
                            solver = solver,
                            multiprocesses_No = multiprocesses_No,
                            method = method,
                            max_times_solve = max_times_solve,
                            verbose= self.verbose   )
        if forPlot:
            self.ResPlotter.add_mz_ranges_to_results(self.res)


    # TODO is the thing below necessary?
    def results_iter(self):
        '''Iterate over results.

        Mainly useful for ETDetective.'''
        for r in self.res:
            for N, info in r["small_graph"].nodes(data=True):
                if info["type"] == "G":
                    mz, intensity, estimate = info['mz'], info['intensity'], info['estimate']
                    L, R = mz.begin, mz.end
                    yield {'L':L,'R':R,'I':intensity,'E':estimate }
        for mz, intensity in izip(*self.spectra['trimmed']):
            prec = self.mz_prec
            yield {'L':mz-prec,'R':mz+prec,'I':intensity,'E':.0 }


    def summarize_results(self):
        '''Summarize the results of MassTodon.'''
        return summarize_results(   peakPicker_stats    = self.peakPicker.stats,
                                    spectra             = self.spectra,
                                    raw_masstodon_res   = self.res              )


    def export_information_for_spectrum_plotting(self, full_info=False):
        '''Provide a generator of dictionaries easy to export as csv file to read in R.'''
        prec = self.mz_prec
        for res in self.ResPlotter.G_info_iter(full_info):
            yield res
        for mz_int in zip(*self.spectra['original']):
            if not mz_int in self.peakPicker.Used_Exps:
                mz, intensity = mz_int
                yield { 'mz_L': mz - prec,
                        'mz_R': mz + prec,
                        'tot_estimate': 0.0,
                        'tot_intensity':intensity,
                        'where': 'not_explainable' }


    # TODO is the thing below necessary?
    def flatten_results(self, minimal_estimated_intensity=100.0):
        '''Return one list of results, one list of difficult cases, and the error.'''
        optimal     = []
        nonoptimal  = []
        total_error  = 0.0
        for mols, error, status in self.res:
            if status=='optimal':
                total_error += error
                for mol in mols:
                    if mol['estimate'] > minimal_estimated_intensity:
                        mol_res = {}
                        for key in ['estimate', 'molType', 'q', 'g', 'formula']:
                            mol_res[key] = mol[key]
                        optimal.append(mol_res)
            else:
                nonoptimal.append(mols)
        return optimal, nonoptimal, total_error


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


    def analyze_reactions(  self,
                            analyzer = 'intermediate',
                            accept_nonOptimalDeconv = False,
                            min_acceptEstimIntensity = 0.0, # might change it to self.spectra['cut_off']
                            verbose=False,
                            **advanced_args ):
        '''Estimate reaction constants and quantities of fragments.'''

        chosen_analyzer = {
            'basic':    analyzer_basic,
            'intermediate':    analyzer_intermediate,
            'advanced': analyzer_advanced
        }[analyzer](self.res, self.Q, self.fasta,
                    accept_nonOptimalDeconv,
                    min_acceptEstimIntensity, verbose )

        return chosen_analyzer.pair()




def MassTodonize(
        fasta,
        precursor_charge,
        mz_prec,
        cut_off         = None,
        opt_P           = None,
        spectrum        = None,
        spectrum_path   = None,
        modifications   = {},
        bootstrap        = False,
        frag_type       = 'cz',
        joint_probability_of_envelope   = .999,
        min_prob_of_envelope_in_picking = .7,
        iso_masses       = None,
        iso_probs        = None,
        L1_x = .001,
        L2_x = .001,
        L1_alpha= .001,
        L2_alpha= .001,
        solver  = 'sequential',
        multiprocesses_No = None,
        method  = 'MSE',
        max_times_solve = 10,
        forPlot = False,
        raw_data= False,
        verbose = False
    ):
    '''Run a full session of MassTodon on your problem.'''
    M = MassTodon(  fasta,
                    precursor_charge,
                    mz_prec,
                    modifications,
                    frag_type,
                    joint_probability_of_envelope,
                    iso_masses,
                    iso_probs,
                    verbose )

    M.read_n_preprocess_spectrum(   spectrum_path,
                                    spectrum,
                                    cut_off,
                                    opt_P   )

    M.run(  solver  = solver,
            multiprocesses_No = multiprocesses_No,
            method  = method,
            max_times_solve = max_times_solve,
            min_prob_per_molecule = min_prob_of_envelope_in_picking,
            bootstrap    = bootstrap,
            forPlot     = forPlot,
            L1_x        = L1_x,
            L2_x        = L2_x,
            L1_alpha    = L1_alpha,
            L2_alpha    = L2_alpha,
            verbose     = verbose       )

    Results = {}
    Results['summary']              = M.summarize_results()
    Results['basic analysis']       = M.analyze_reactions('basic')
    Results['intermediate analysis']= M.analyze_reactions('intermediate')
    Results['advanced analysis']    = M.analyze_reactions('advanced')
    if raw_data:
        Results['raw estimates'] = M.res
    if forPlot:
        Results['short data to plot']   = M.export_information_for_spectrum_plotting(False)
        Results['long data to plot']    = M.export_information_for_spectrum_plotting(True)
        Results['original spectrum']    = M.spectrum_iter('original')
    return Results
