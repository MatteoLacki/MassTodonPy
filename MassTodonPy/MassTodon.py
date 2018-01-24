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
from __future__ import absolute_import, division, print_function
from collections import Counter
from math import ceil, log10
import csv

from MassTodonPy.Data.Constants import eps
from MassTodonPy.MatchMaker.CzMatch import CzMatch
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch
from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Precursor.Precursor import Precursor
from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.Deconvolution.Deconvolve import deconvolve


class MassTodon(object):
    def __init__(self,
                 spectrum,
                 precursor,
                 mz_precision,
                 preprocessing_args={},
                 isospec_args={},
                 deconvolution_args={},
                 solver_args={},
                 simple_cz_match=False,
                 _devel=False):
        """Run a full session of the MassTodon.

        Parameters
        ==========
        spectrum : string, tuple of numpy arrays, or Spectrum
            The string path can end with:
                *.txt, for spectra saved in tab separated format.
                *.mzXml, for spectra saved in mzXml format.
                *.mzml, for spectra saved mzml format.
            The tuple consists of two numpy arrays:
                one with m over z ratios,
                the other with intensities.
            Details of the Spectrum can be found here [TODO add link].
        precursor : dict
            A dictionary with fasta and charge. Also accepts the precursor's name,
            modifications, fragmentation_type, blocked_fragments, and distance_charges.
            Details can be found here.
        mz_precision : float
            The precision in the m/z domain: how far away [in Daltons] should
            one search for peaks from the theoretical mass.
        preprocessing_args : dict
            Arguments for spectrum preprocessing, such as min_prob_per_molecule and mz_precision.
        deconvolution_args : dict
            Arguments for the deconvolution stage, such as: method, mz_tol, min_prob_per_molecule
        solver_args : dict
            Arguments for the solver.
        """

        assert isinstance(mz_precision, float) and mz_precision >= 0.0
        # precision one digit lower than the input precision of spectra, eg.
        # mz_precision = .05  --> prec_digits = 3
        # mz_precision = .005 --> prec_digits = 4
        mz_digits_tmp = int(ceil(-log10(mz_precision)))
        self.mz_digits = deconvolution_args.get('mz_digits', mz_digits_tmp)
        self.precursor = Precursor(**precursor)
        self.spectrum = Spectrum(spectrum=spectrum,
                                 mz_digits=self.mz_digits,
                                 **preprocessing_args)

        isospec_args['mz_digits'] = self.mz_digits
        self.minimal_intensity = preprocessing_args.get('minimal_intensity', eps)
        self._solutions = deconvolve(self.precursor.molecules(),
                                     self.spectrum,
                                     isospec_args=isospec_args,
                                     solver_args=solver_args,
                                     **deconvolution_args)
        if _devel:
            self._solutions = list(self._solutions)
        self._raw_estimates = list(self.get_raw_estimates(minimal_intensity=eps)) #TODO terminate
        if simple_cz_match:
            #TODO adjust input reading
            self.simple_cz_match = SimpleCzMatch(self._raw_estimates, self.precursor)
        self.cz_match = CzMatch(self._raw_estimates, self.precursor)

    def get_raw_estimates(self, minimal_intensity=eps):
        """Iterate over estimates with intensity greater than the minimal_intensity."""
        for sol in self._solutions:
            res = sol.report()
            if res['status'] is not 'ValueError':
                for mol in res['alphas']:
                    estimate = int(mol['estimate'])
                    if estimate >= minimal_intensity:
                        mol = mol['molecule']
                        yield mol, estimate

    def get_estimates_for_plot(self, minimal_intensity=eps):
        """Probably just like function above."""
        pass

    def write(self, path):
        """Write the spectrum to a csv or tsv file.

        Parameters
        ==========
        path : str
            A path to the file to write to.
        """
        file_path, file_name, file_ext = parse_path(path)
        assert file_ext in ('.csv', '.tsv'), "Writing only to csv or tsv."
        delimiter = ',' if file_ext == '.csv' else '\t'

        # write intensities of molecules
        path_intensities = path.replace('.', '_intensities.')
        with open(path_intensities, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(['molecule', 'source', 'formula', 'charge', 'quenched charge',
                             'monoisotopic m/z', 'average m/z', 'intensity'])
            for m, intensity in self.get_estimates():
                writer.writerow([m.name, m.source, str(m.formula), m.q, m.g,
                                 m.monoisotopic_mz, m.mean_mz, intensity])

        # write intensities of reactions
        path_reactions = path.replace('.', '_reactions.')
        with open(path_reactions, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(['reaction', 'intensity'])
            for name, intensity in self.cz_match.intensities.items():
                if isinstance(intensity, (Counter, dict)):
                    for k, v in intensity.items():
                        writer.writerow([name + str(k), v])
                else:
                    writer.writerow([name, intensity])

        # write probabilities of reactions
        path_probabilities = path.replace('.', '_probabilities.')
        with open(path_probabilities, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(['reaction', 'probability'])
            for name, intensity in self.cz_match.probabilities.items():
                if isinstance(intensity, (Counter, dict)):
                    for k, v in intensity.items():
                        writer.writerow([name + '_' + str(k), v])
                else:
                    writer.writerow([name, intensity])

    def plot(self, simple=True):
        pass
