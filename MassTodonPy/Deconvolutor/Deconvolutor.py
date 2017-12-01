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
from __future__ import absolute_import, division, print_function
# from multiprocessing import Pool
from itertools import repeat
from six.moves import zip
# import traceback  # to see traceback for children processes

from MassTodonPy.Deconvolutor.DeconvolutionProblem import DeconvolutionProblem
from MassTodonPy.Deconvolutor.PeakPicker import get_deconvolution_graphs

# TODO play with multiprocessing later on.

# def worker(problem):
    # graph, deconv_args = args
    # problem = DeconvolutionProblem(graph, **deconv_args)
    # problem.solve()
    # return problem

def deconvolve(molecules,
               spectrum,
               # processes=0,
               **deconvolution_args):
    """Solve deconvolution problems."""
    for graph in get_deconvolution_graphs(molecules, spectrum):
        problem = DeconvolutionProblem(graph, **deconvolution_args)
        problem.solve()
        yield problem

    # if processes in (0, 1):
    #     out = [worker((g, deconvolution_args)) for g in graphs]
    # else:
    #     assert processes > 1, "processes in {0, 1, 2, ...}"
    #     graphs.sort(reverse=True, key=len)
    #     P = Pool(int(processes))
    #     out = P.map(worker,
    #                 zip(graphs,
    #                     repeat(deconvolution_args)))
    #     P.close()
    #     P.join()
    # return problems
