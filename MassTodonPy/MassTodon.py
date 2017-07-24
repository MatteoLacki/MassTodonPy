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


# Here is our little pet MassTodon!
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

#### MassTodonPy modules
from Formulator         import make_formulas
from IsotopeCalculator  import IsotopeCalculator
from PeakPicker         import PeakPicker
from Solver             import solve
from Parsers            import read_n_preprocess_spectrum
from MatchMaker         import match_cz_ions
from Visualization      import ResultsPlotter, make_highcharts
from Summarator         import summarize_results
from Outputing          import write_raw_to_csv, write_counts_n_probs_to_csv, write_summary_to_csv

#### standard modules
from itertools          import izip, chain
from math               import ceil, log10
from intervaltree       import Interval as interval, IntervalTree
from time               import time


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
        self.small_graphs_no_G = None

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
            for_plot            = False,
            **args ):
        '''Perform the deconvolution of problems.'''

        self.problems = self.peakPicker.get_problems(
            massSpectrum            = self.spectra['trimmed'],
            min_prob_per_molecule   = min_prob_per_molecule )

        self.spectra['intensity of peaks paired with isotopologues'] = self.peakPicker.stats['total intensity of experimental peaks paired with isotopologues']

        self.res, self.solver_stats = solve(
                            problems = self.problems,
                            args   = args,
                            solver = solver,
                            multiprocesses_No = multiprocesses_No,
                            method = method,
                            max_times_solve = max_times_solve,
                            verbose= self.verbose   )

        if self.verbose:
            print 'Solver stats:'
            print self.solver_stats
            print

        if for_plot:
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
        return summarize_results(   spectra             = self.spectra,
                                    raw_masstodon_res   = self.res              )


    #TODO use Bokeh from Python directly!
    def make_data_for_spectrum_plot(self):
        '''Provide estimates easily exportable to Bokeh via R.'''
        data_4_plot = self.ResPlotter.make_data_for_spectrum_plot()
        data_4_plot['remaining_peaks'] = [
            {   'mz_L':         mz - self.mz_prec,
                'mz_R':         mz + self.mz_prec,
                'tot_estimate': 0.0,
                'tot_intensity':intensity
            }   for mz, intensity in izip(*self.spectra['original'])
                if not (mz, intensity) in self.peakPicker.Used_Exps
            ]


        return data_4_plot


    def analyze_reactions(  self,
                            analyzer                = 'intermediate',
                            accept_nonOptimalDeconv = False,
                            min_acceptEstimIntensity= 0.0, # might change it to self.spectra['cut_off']
                            verbose                 = False,
                            **advanced_args     ):

        '''Estimate reaction constants and quantities of fragments.'''
        return match_cz_ions(   results_to_pair         = self.res,
                                Q                       = self.Q,
                                fasta                   = self.fasta,
                                min_acceptEstimIntensity= min_acceptEstimIntensity,
                                analyzer                = analyzer,
                                accept_nonOptimalDeconv = accept_nonOptimalDeconv,
                                verbose                 = verbose,
                                advanced_args           = advanced_args   )



def MassTodonize(
        fasta,
        precursor_charge,
        mz_prec,
        cut_off         = None,
        opt_P           = None,
        spectrum        = None,
        spectrum_path   = None,
        modifications   = {},
        frag_type       = 'cz',
        joint_probability_of_envelope   = .999,
        min_prob_of_envelope_in_picking = .7,
        iso_masses  = None,
        iso_probs   = None,
        L1_x = .001,
        L2_x = .001,
        L1_alpha= .001,
        L2_alpha= .001,
        solver  = 'sequential',
        multiprocesses_No = None,
        method  = 'MSE',
        max_times_solve = 10,
        for_plot= False,
        highcharts = False,
        raw_data= False,
        analyze_raw_data = True,
        output_csv_path  = None,
        output_deconvolution_threshold = 0.0,
        etdetective = False,
        verbose = False
    ):
    '''Run a full session of MassTodon on your problem.'''

    T0 = time()
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

    M.run(  solver = solver,
            multiprocesses_No = multiprocesses_No,
            method = method,
            max_times_solve = max_times_solve,
            min_prob_per_molecule = min_prob_of_envelope_in_picking,
            for_plot  = for_plot,
            L1_x      = L1_x,
            L2_x      = L2_x,
            L1_alpha  = L1_alpha,
            L2_alpha  = L2_alpha,
            verbose   = verbose    )

    results = {}

    if analyze_raw_data:
        results['summary']              = M.summarize_results()
        results['basic_analysis']       = M.analyze_reactions('basic')
        results['intermediate_analysis']= M.analyze_reactions('intermediate')
        results['advanced_analysis']    = M.analyze_reactions('advanced')

    if verbose:
        print 'L1_error_value_error/intensity_within_tolerance', results['summary']['L1_error_value_error/intensity_within_tolerance']
        print

    if raw_data:
        results['raw_estimates'] = M.res
        results['spectra'] = M.spectra

    if for_plot:
        results['for_plot']= M.make_data_for_spectrum_plot()

    if highcharts:
        algos = {   'basic_analysis':
        results['basic_analysis'],
                    'intermediate_analysis':    results['intermediate_analysis'],
                    'advanced_analysis':        results['advanced_analysis']        }

        results['highcharts'] = make_highcharts(
            fasta   = fasta,
            Q       = precursor_charge,
            raw_estimates = M.res,
            algos   = algos)

    T1 = time()
    results['summary']['total_time'] = T1-T0

    if output_csv_path:
        write_raw_to_csv(M.res, output_csv_path, output_deconvolution_threshold)
        write_counts_n_probs_to_csv(results, fasta, output_csv_path)
        write_summary_to_csv(results, output_csv_path)

    if etdetective:
        results['etdetective'] = results_to_etdetective( M.res, M.fasta, M.modifications )

    if verbose:
        print 'Total analysis took', T1-T0

    return results
