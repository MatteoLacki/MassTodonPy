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

from MassTodonPy.Misc.strings import add_backslash

def add_spectra_plots_parsing(parser):
    """Add parsing of spectra plots."""
    parser.add_argument("--no_spectrum_plot",
                        dest='spectrum_plot',
                        action='store_const',
                        const=False,
                        default=True,
                        help="Skip the spectrum plot.")
    parser.add_argument("--show_spectrum_plot",
                        dest='show_spectrum_plot',
                        action='store_const',
                        const=True,
                        default=False,
                        help="Show spectrum plot.")
    parser.add_argument("-spectrum_plot_width",
                        dest='width_spectrum_plot',
                        default=800,
                        type=int,
                        help="Width of the spectrum plot.")
    parser.add_argument("-spectrum_plot_height",
                        dest='height_spectrum_plot',
                        default=600,
                        type=int,
                        help="Height of the spectrum plot.")
    # aggregated precursors plot
    parser.add_argument("--no_aggregated_precursors_plot",
                        dest='aggregated_precursors_plot',
                        action='store_const',
                        const=False,
                        default=True,
                        help="Skip the plot of aggregated precursors.")
    parser.add_argument("--show_aggregated_precursors_plot",
                        dest='show_aggregated_precursors_plot',
                        action='store_const',
                        const=True,
                        default=False,
                        help="Do not show the aggregated precursors plot.")
    parser.add_argument("-aggregated_precursors_plot_width",
                        dest='width_aggregated_precursors_plot',
                        default=800,
                        type=int,
                        help="Width of the plot of aggregated precursors.")
    parser.add_argument("-spectrum_height",
                        dest='height_aggregated_precursors_plot',
                        default=600,
                        type=int,
                        help="Height of the plot of aggregated precursors.")
    # aggregated fragments plot
    parser.add_argument("--no_aggregated_fragments_plot",
                        dest='aggregated_fragments_plot',
                        action='store_const',
                        const=False,
                        default=True,
                        help="Skip the aggregated fragments plot.")
    parser.add_argument("--show_aggregated_fragments_plot",
                        dest='show_aggregated_fragments_plot',
                        action='store_const',
                        const=True,
                        default=False,
                        help="Do not show the aggregated fragments plot.")
    parser.add_argument("-aggregated_fragments_plot_width",
                        dest='width_aggregated_fragments_plot',
                        help="Width of the aggregated fragments plot.",
                        default=0,
                        type=int)
    parser.add_argument("-aggregated_fragments_plot_height",
                        dest='height_aggregated_fragments_plot',
                        help="Height of the aggregated fragments plot.",
                        default=0,
                        type=int)
    # estimated aggregated fragments plot
    parser.add_argument("--no_estimated_aggrageted_fragments_plot",
                        dest='estimated_aggregated_fragments_plot',
                        action='store_const',
                        const=False,
                        default=True,
                        help="Skip the estimated aggregated fragments plot.")
    parser.add_argument("--show_estimated_aggregated_fragments_plot",
                        dest='show_estimated_aggregated_fragments_plot',
                        action='store_const',
                        const=True,
                        default=False,
                        help="Do not show the aggregated fragments plot.")
    parser.add_argument("-estimated_aggregated_fragments_plot_width",
                        dest='width_estimated_aggrageted_fragments_plot',
                        help="Width of the estimated aggregated fragments plot.",
                        default=0,
                        type=int)
    parser.add_argument("-estimated_aggregated_fragments_plot_height",
                        dest='height_estimated_aggregated_fragments_plot',
                        help="Height of the estimated aggregated fragments plot.",
                        default=0,
                        type=int)


def add_verbosity_parsing(parser):
    parser.add_argument("--verbose",
                        dest='verbose',
                        action='store_const',
                        const=True,
                        default=False,
                        help="Print out the messages from MassTodon.")

def add_max_times_parsing(parser):
    parser.add_argument("-max_times",
                        dest="_max_times",
                        help="The maximal number of times to run CVXOPT.",
                        type=int,
                        default=10)

def add_output_parsing(parser):
    parser.add_argument("-o",
                        "--output_path",
                        dest='output',
                        help="Path to the output.\
                              If not provided, spectrum is used path by default.",
                        default='output/',
                        type=add_backslash)
