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
from __future__ import absolute_import

from MassTodonPy.Deconvolutor.DeconvolutionProblem import DeconvolutionProblem
from MassTodonPy.Deconvolutor.PeakPicker import get_deconvolution_problems as get_graphs


def deconvolve(molecules,
               spectrum,
               method="Matteo",
               isospec_args={},
               mz_tol_args={},
               deconvolution_args={},
               **args):
    graphs = get_graphs(molecules,
                        spectrum, 
                        method,
                        isospec_args,
                        mz_tol_args)
    """Solve the deconvolution problems."""
    for graph in graphs:
        if method == 'Matteo':
            problem = DeconvolutionProblem(graph, **deconvolution_args)
            problem.solve()
            yield problem
        else:
            yield graph
        