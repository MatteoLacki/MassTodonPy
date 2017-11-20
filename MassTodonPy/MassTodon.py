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

# from six.moves import zip
from math import ceil, log10
# from time import time

from MassTodonPy.MoleculeMaker.Precursor import Precursor
from MassTodonPy.MoleculeMaker.MoleculeMaker import get_molecules
from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.Data.Constants import eps


# from MassTodonPy.MoleculeMaker.MoleculeMaker import make_molecules
# from IsotopeCalculator import IsotopeCalculator

# from PeakPicker.peakPicker import PeakPicker
# from MatchMaker.match_cz_ions import match_cz_ions
# from Outputting.export_outputs import OutputExporter
# from Outputting.write_to_csv import write_raw_to_csv,
# write_counts_n_probs_to_csv, write_summary_to_csv
# from Solver.solver import solve
# from Summarator.summarator import summarize_results


def MassTodonize(precursor_name,
                 precursor_fasta,
                 precursor_charge,
                 precursor_modifications={},
                 m_over_z_precision=0.05,
                 fragmentation_type='cz',
                 blockedFragments=set(['c0']),
                 minimal_distance_between_charges=5,
                 joint_probability_of_the_envelope=0.999,
                 min_prob_of_envelope_in_picking=.7,
                 spectrum='',
                 spectrum_minimal_intensity=eps,
                 spectrum_percent_top_peaks=1.0,
                 processes_no=1,
                 _faster_mzxml=False,
                 _max_times_solve=10,
                 _isotope_masses=None,
                 _isotope_probabilities=None,
                 _fitting_criterion='MSE',
                 _L1_x=.001,
                 _L2_x=.001,
                 _L1_alpha=.001,
                 _L2_alpha=.001,
                 _verbose=False):
    """
    Run a full session of MassTodon on your problem.
    Parameters
    ==========
    precursor_name : str
        The name of the input substance.

    precursor_fasta : str
        The fasta of the input substance.

    precursor_charge : int
        The charge of the input substance.

    precursor_modifications :
        The key in the modifications' dictionary should with of form
        **<N|Calpha|C><amino acid No>**, e.g. C10.
        The values contain dictionaries with diffs in elements,
        e.g. **{'H': 1, 'N': 1, 'O': -1}**. TODO: check it.

    m_over_z_precision : float
        The precision in the m/z domain: how far away [in Daltons] should
        one search for peaks from the theoretical mass.

    fragmentation_type : str
        Types of fragments you want to find.
        Warning
        =======
        Supports c and z fragmentations. So you can input "cz", "c", or "z".

    blockedFragments : set of strings
        Which fragments should not appear. Defaults to 'c0'.
        MoleculeMaker adds to this all proline-blocked fragments.

    minimal_distance_between_charges: int
        The minimal distance between charges on the protein.
        If set to 5, at least 4 amino acids must la
        between consecutive *charged* amino acids.

    joint_probability_of_envelope : float
        The joint probability threshold for generating
        theoretical isotopic envelopes.

    min_prob_of_envelope_in_picking : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.

    spectrum : path string or a tuple of numpy arrays (masses, intensities).

    spectrum_minimal_intensity : float
        Experimental peaks with lower height will be trimmed.

    spectrum_percent_top_peaks : float
        The percentage of the heighest experimental peaks used in the analysis.

    processes_no : int
        How much processes should we use to perform the fitting?
        One process is the same as sequential calculations.

    _max_times_solve : int
        How many times the CVXOPT solver should try to find the optimal
        deconvolution of the isotopic envelopes.

        Note
        ====
        CVXOPT sometimes does not meet its optimality criteria.
        However, it uses underneath the OPENBLAS library and this baby has
        its own dynamic scheduling that might sometimes influence the output
        of the calculations.
        Here, we try to use it to our advantage, in order to help CVXOPT
        to meet its optimisation criteria.

        This is clearly not the nice way to do it, but we will soon go bayesian
        anyway and use a different approach to perform the deconvolution.

    _isotope_masses : dict
        The isotopic masses, e.g. from IUPAC.

    _isotope_probabilities : dict
        The isotopic frequencies, e.g. from IUPAC.

    _fitting_criterion : str
        The minimization criterion used.
        Warning
        =======
            For now, only mean square error (MSE) considered.

    _L1_x : float
        The $L_1$ penalty for flows (intensities attributed to particular
        isotopologous and experimental groups) in the deconvolution problem.

    _L2_x : float
        The $L_2$ penalty for flows (intensities attributed to particular
        isotopologous and experimental groups) in the deconvolution problem.

    _L1_alpha : float
        The $L_1$ penalty for total intensities attributed to particular
        molecular species in the deconvolution problem.

    _L2_alpha : float
        The $L_2$ penalty for total intensities attributed to particular
        molecular species in the deconvolution problem.

    _verbose : boolean
        Do you want to check what MassTodon did undercover?

    Returns
    =======
    out : a results object
    """
    assert isinstance(m_over_z_precision, float) and m_over_z_precision >= 0.0

    # precision one digit lower than the input precision of spectra, eg.
    # m_over_z_precision = .05     -->     prec_digits = 3
    # m_over_z_precision = .005    -->     prec_digits = 4
    mz_precision_digits = int(ceil(-log10(m_over_z_precision)))

    precursor = Precursor(name=precursor_name,
                          fasta=precursor_fasta,
                          q=precursor_charge,
                          modifications=precursor_modifications)

    molecules = get_molecules([precursor],
                              blockedFragments,
                              fragmentation_type,
                              minimal_distance_between_charges)

    spectrum = Spectrum(spectrum,
                        mz_precision_digits,
                        spectrum_minimal_intensity,
                        spectrum_percent_top_peaks,
                        _faster_mzxml)
    pass

#
#
# class MassTodon():
#     """Make MassTodons less extinct by giving to the World its instance.
#
#     Parameters
#     ----------
#     fasta : str
#         The fasta of the studied molecular species.
#
#     precursor_charge : int
#         The charge of the precursor ion.
#
#     m_over_z_precision: float
#         The precision in the m/z domain: how far away [in Daltons]
#         should we search for peaks from the theoretical mass.
#
#         For instance, consider H20 and m_over_z_precision = 0.05 Da.
#         Then, the mass of the lightest isotopic variant equals
#         1.0078 x 2 + 15.9949 = 18.0105 Da.
#         Then, intensity within a range 18.0105 ± 0.05 Da will be considered
#         to potentially originate from the water molecules.
#
#     modifications : dict
#         A dictionary of modifications.
#
#         The key in the modifications' dictionary should with of form
#         **<N|Calpha|C><amino acid No>**, e.g. C10.
#
#         The values contain dictionaries with diffs in elements,
#         e.g. **{'H': 1, 'N': 1, 'O': -1}**.
#
#     minimal_distance_between_charges : int
#         The minimal distance between charges on the protein.
#         If set to 5, at least 4 amino acids must lay between
#         consecutive *charged* amino acids.
#     """
#
#     def __init__(self,
#                  precursor_name,
#                  precursor_fasta,
#                  precursor_charge,
#                  precursor_modifications={},
#                  m_over_z_precision=0.05,
#                  fragmentation_type='cz',
#                  blockedFragments=set(['c0']),
#                  minimal_distance_between_charges=5,
#                  joint_probability_of_the_envelope=0.999,
#                  _isotope_masses=None,
#                  _isotope_probabilities=None,
#                  _amino_acids=None,
#                  _verbose=False):
#         self.mz_prec = m_over_z_precision
#         # precision one digit lower than the input precision of spectra, eg.
#         # mz_prec = .05     -->     prec_digits = 3
#         # mz_prec = .005    -->     prec_digits = 4
#         self.prec_digits = int(ceil(-log10(self.mz_prec)))
#
#         self.precursor = Precursor(name=precursor_name,
#                                    fasta=precursor_fasta,
#                                    q=precursor_charge,
#                                    modifications=precursor_modifications)
#
#         self.molecules = molecules([self.precursor],
#                                    blockedFragments,
#                                    fragmentation_type,
#                                    minimal_distance_between_charges)
#
#         self.spectra = None
#         self.verbose = _verbose
#
#     def read_n_preprocess_spectrum(self,
#                                    path=None,
#                                    spectrum=None,
#                                    cut_off=None,
#                                    opt_P=None):
#         pass

    # #  TODO is the thing below necessary?
    # def spectrum_iter(self, spectrum_type):
    #     assert spectrum_type in ['original', 'trimmed'], "No such kind of spectrum: %s." % spectrum_type
    #     for mz, intensity in zip(*self.spectra[spectrum_type]):
    #         yield {'mz': mz, 'intensity': intensity}

    # def run(self,
    #         solver='sequential',
    #         multiprocesses_No=None,
    #         method='MSE',
    #         max_times_solve=5,
    #         min_prob_per_molecule=.75,
    #         for_plot=False,
    #         **args):
    #     '''Perform the deconvolution of problems.'''
    #
    #     self.spectra = read_n_preprocess_spectrum(path=path,
    #                                               spectrum=spectrum,
    #                                               prec_digits=self.prec_digits,
    #                                               cut_off=cut_off,
    #                                               opt_P=opt_P,
    #                                               verbose=self.verbose)

        # self.molecules = make_molecules(fasta=self.fasta,
        #                              Q=self.Q,
        #                              amino_acids=self.amino_acids,
        #                              distance_charges=self.distance_charges,
        #                              modifications=self.modifications)
        #
        # self.IsoCalc = IsotopeCalculator(jP=joint_probability_of_envelope,
        #                                  prec_digits=self.prec_digits,
        #                                  iso_masses=iso_masses,
        #                                  iso_probs=iso_probs,
        #                                  verbose=verbose)
        #
        # self.peakPicker = PeakPicker(
        #     _Forms=self.formulas, _IsoCalc=self.IsoCalc,
        #     mz_prec=mz_prec, verbose=verbose)
        #
        # self.ResPlotter = OutputExporter(mz_prec)
        #
        # self.problems = self.peakPicker.get_problems(
        #     spectrum=self.spectra['trimmed'],
        #     min_prob_per_molecule=min_prob_per_molecule)
        #
        # self.spectra['intensity of peaks paired with isotopologues'] = self.peakPicker.stats['total intensity of experimental peaks paired with isotopologues']
        #
        # self.res, self.solver_stats = solve(
        #                     problems=self.problems,
        #                     args=args,
        #                     solver=solver,
        #                     multiprocesses_No=multiprocesses_No,
        #                     method=method,
        #                     max_times_solve=max_times_solve,
        #                     verbose=self.verbose)
        #
        # if self.verbose:
        #     print('Solver stats:')
        #     print(self.solver_stats)
        #     print()
        #
        # if for_plot:
        #     self.ResPlotter.add_mz_ranges_to_results(self.res)

    # # TODO is the thing below necessary?
    # def results_iter(self):
    #     '''Iterate over results.
    #     Mainly useful for ETDetective.'''
    #     for r in self.res:
    #         for N, info in r["small_graph"].nodes(data=True):
    #             if info["type"] == "G":
    #                 mz, intensity, estimate = info['mz'], info['intensity'], info['estimate']
    #                 L, R = mz.begin, mz.end
    #                 yield {'L': L, 'R': R, 'I': intensity, 'E': estimate}
    #     for mz, intensity in zip(*self.spectra['trimmed']):
    #         prec = self.mz_prec
    #         yield {'L': mz-prec, 'R': mz+prec, 'I': intensity, 'E': .0}
    #
    # def summarize_results(self):
    #     '''Summarize the results of MassTodon.'''
    #     return summarize_results(spectra=self.spectra,
    #                              raw_masstodon_res=self.res)
    #
    # # TODO use Bokeh from Python directly!
    # def make_data_for_spectrum_plot(self):
    #     '''Provide estimates easily exportable to Bokeh via R.'''
    #     data_4_plot = self.ResPlotter.make_data_for_spectrum_plot()
    #     data_4_plot['remaining_peaks'] = [
    #         {'mz_L': mz - self.mz_prec,
    #          'mz_R': mz + self.mz_prec,
    #          'tot_estimate': 0.0,
    #          'tot_intensity': intensity}
    #         for mz, intensity in zip(*self.spectra['original'])
    #         if not (mz, intensity) in self.peakPicker.Used_Exps]
    #     return data_4_plot
    #
    # def analyze_reactions(self,
    #                       analyzer='intermediate',
    #                       min_acceptEstimIntensity=0.0,  # self.spectra['cut_off']
    #                       verbose=False,
    #                       **advanced_args):
    #     '''Estimate reaction constants and quantities of fragments.'''
    #     return match_cz_ions(results_to_pair=self.res,
    #                          Q=self.Q,
    #                          fasta=self.fasta,
    #                          min_acceptEstimIntensity=min_acceptEstimIntensity,
    #                          analyzer=analyzer,
    #                          verbose=verbose,
    #                          advanced_args=advanced_args)

#
# def MassTodonize(
#         fasta,
#         precursor_charge,
#         mz_prec,
#         cut_off=None,
#         opt_P=None,
#         spectrum=None,
#         spectrum_path=None,
#         modifications={},
#         frag_type='cz',
#         distance_charges=5,
#         joint_probability_of_envelope=.999,
#         min_prob_of_envelope_in_picking=.7,
#         iso_masses=None,
#         iso_probs=None,
#         L1_x=.001,
#         L2_x=.001,
#         L1_alpha=.001,
#         L2_alpha=.001,
#         solver='sequential',
#         multiprocesses_No=None,
#         method='MSE',
#         max_times_solve=10,
#         for_plot=False,
#         raw_data=False,
#         analyze_raw_data=True,
#         output_csv_path=None,
#         output_deconvolution_threshold=0.0,
#         verbose=False):
#     """Run a full session of MassTodon on your problem.
#
#     Parameters
#     ==========
#     fasta : str
#         The fasta of the studied molecular species.
#
#     precursor_charge : int
#         The charge of the precursor ion.
#
#     mz_prec : float
#         The precision in the m/z domain: how far away [in Daltons] should
#         one search for peaks from the theoretical mass.
#
#     cut_off : float
#         The cut off value for peak intensity.
#
#         Warning
#         =======
#         You should provide either the path or the spectrum,
#         not both simultaneously.
#
#     opt_P :
#         The percentage of the heighest peaks being used in the analysis.
#
#         Warning
#         =======
#         You should provide either the path or the spectrum,
#         not both simultaneously.
#
#     spectrum : tuple of numpy arrays
#         A mass spectrum that you prepared *30 minutes before the programme*,
#         consisting of a tuple of numpy arrays.
#         The first array corresponds to different **m/z** ratios,
#         the other contains their respective intensity values.
#
#         Warning
#         =======
#         You should set either *spectrum* or *spectrum_path*,
#         but not both simultaneously.
#
#     spectrum_path : str
#         The path to spectrum in either mzXml or tsv format.
#
#         Warning
#         =======
#         You should set either *spectrum* or *spectrum_path*,
#         but not both simultaneously.
#
#     modifications : dict
#         A dictionary of modifications.
#
#         The key in the modifications' dictionary should with of form
#         **<N|Calpha|C><amino acid No>**, e.g. C10.
#         The values contain dictionaries with diffs in elements,
#         e.g. **{'H': 1, 'N': 1, 'O': -1}**.
#
#     frag_type : str
#         The type of framgentation.
#
#         Warning
#         =======
#         Only **cz** fragmentation supported.
#
#     distance_charges : int
#         The minimal distance between charges on the protein.
#         If set to 5, at least 4 amino acids must la
#         between consecutive *charged* amino acids.
#
#     joint_probability_of_envelope : float
#         The joint probability threshold for generating
#         theoretical isotopic envelopes.
#
#     min_prob_of_envelope_in_picking : float
#         The minimal probability an envelope has to scoop
#         to be included in the deconvolution graph.
#
#     iso_masses : dict
#         The isotopic masses.
#
#     iso_probs : dict
#         The isotopic frequencies.
#
#     L1_x : float
#         The $L_1$ penalty for flows (intensities attributed to particular
#         isotopologous and experimental groups) in the deconvolution problem.
#
#     L2_x : float
#         The $L_2$ penalty for flows (intensities attributed to particular
#         isotopologous and experimental groups) in the deconvolution problem.
#
#     L1_alpha : float
#         The $L_1$ penalty for total intensities attributed to particular
#         molecular species in the deconvolution problem.
#
#     L2_alpha : float
#         The $L_2$ penalty for total intensities attributed to particular
#         molecular species in the deconvolution problem.
#
#     solver : str
#         Select the mode of solving the deconvolution tasks:
#         either *sequential* or *multiprocessing*.
#
#         Note
#         ====
#         Despite the possibility to solve simultaneously deconvolution problems
#         for different regions of the mass spectrum, the sequential solver is
#         usually faster due to overheads on process comunication or
#         my terrible coding skills.
#
#         Warning
#         =======
#         While implementing some multiprocessing routine of your own that could
#         involve MassTodonPy, please use rather the sequential fitting.
#         Altogether, making CVXOPT and multiprocessing work is a nightmare.
#
#     multiprocesses_No : int
#         How many processes should be used during the deconvolution.
#         Defaults to the number of cores.
#
#     method : str
#         What kind of optimization scheme should be used for deconvolution?
#
#         Warning
#         =======
#         Currently only accepts "MSE", which means that the deconvolution
#         will be performed using the Means Square Error minimization,
#         which leads to a quadratic programme with linear constraints.
#
#     max_times_solve : int
#         How many times the CVXOPT solver should try to find the optimal
#         deconvolution of the isotopic envelopes.
#
#         Note
#         ====
#         CVXOPT sometimes does not meet its optimality criteria.
#         However, it uses underneath the OPENBLAS library and this baby has
#         its own dynamic scheduling that might sometimes influence the output
#         of the calculations.
#         Here, we try to use it to our advantage, in order to help CVXOPT
#         to meet its optimisation criteria.
#
#         This is clearly not the nice way to do it, but we will soon go bayesian
#         anyway and use a different approach to perform the deconvolution.
#
#     for_plot : boolean
#         Add data for plot to results.
#
#     raw_data : boolean
#         Should we add the raw output (peak assignement) to the output?
#
#     analyze_raw_data : boolean
#         Select **True** to pair the estimates of fragments and perform their
#         analysis.
#         Among the output, you will find probabilities of fragmentations and
#         the ETnoD and PTR reactions, as well as their intensities.
#
#     output_csv_path : str
#         Path to where you want to save the output, if you want to save it.
#
#     output_deconvolution_threshold : float
#         Only results with estimates above this threshold will be shown
#         in the final output.
#
#     verbose : boolean
#         Do you want to check what MassTodon did undercover?
#
#
#     Returns
#     =======
#     results : dict
#         A dictionary containing the results of MassTodonPy.
#     """
#
#     T0 = time()
#
#     M = MassTodon(fasta=fasta,
#                   precursor_charge=precursor_charge,
#                   mz_prec=mz_prec,
#                   modifications=modifications,
#                   frag_type=frag_type,
#                   distance_charges=distance_charges,
#                   joint_probability_of_envelope=joint_probability_of_envelope,
#                   iso_masses=iso_masses,
#                   iso_probs=iso_probs,
#                   verbose=verbose)
#
#     M.read_n_preprocess_spectrum(spectrum_path, spectrum,
#                                  cut_off, opt_P)
#
#     M.run(solver=solver,
#           multiprocesses_No=multiprocesses_No,
#           method=method,
#           max_times_solve=max_times_solve,
#           min_prob_per_molecule=min_prob_of_envelope_in_picking,
#           for_plot=for_plot,
#           L1_x=L1_x,
#           L2_x=L2_x,
#           L1_alpha=L1_alpha,
#           L2_alpha=L2_alpha,
#           verbose=verbose)
#
#     results = {}
#
#     if analyze_raw_data:
#         results['summary'] = M.summarize_results()
#         results['basic_analysis'] = M.analyze_reactions('basic')
#         results['intermediate_analysis'] = M.analyze_reactions('intermediate')
#         results['advanced_analysis'] = M.analyze_reactions('advanced')
#
#     if verbose:
#         print('L1_error_value_error/intensity_within_tolerance', results['summary']['L1_error_value_error/intensity_within_tolerance'])
#         print()
#
#     if raw_data:
#         results['raw_estimates'] = M.res
#         results['spectra'] = M.spectra
#
#     if for_plot:
#         results['for_plot'] = M.make_data_for_spectrum_plot()
#
#     # if highcharts:
#     #     algos = {'basic_analysis': results['basic_analysis'],
#     #              'intermediate_analysis': results['intermediate_analysis'],
#     #              'advanced_analysis': results['advanced_analysis']}
#     #
#     #     results['highcharts'] = make_highcharts(fasta=fasta,
#     #                                             Q=precursor_charge,
#     #                                             raw_estimates=M.res,
#     #                                             algos=algos)
#
#     T1 = time()
#     results['summary']['total_time'] = T1-T0
#
#     if output_csv_path:
#         write_raw_to_csv(M.res, output_csv_path, output_deconvolution_threshold)
#         write_counts_n_probs_to_csv(results, fasta, output_csv_path)
#         write_summary_to_csv(results, output_csv_path)
#
#     if verbose:
#         print('Total analysis took', T1-T0)
#
#     return results
