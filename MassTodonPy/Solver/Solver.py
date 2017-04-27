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
from MassTodonPy.Deconvolutor import deconvolve

class Solver(object):
    def __init__(self, problemsGenerator):
        self.prob_gen = problemsGenerator
    def run(self, args, method='MSE'):
        raise NotImplementedError


class SequentialSolver(Solver):
    def run(self, args, method='MSE'):
        results = []
        for SFG in self.prob_gen:
            res = deconvolve(   SFG     = SFG,
                                args    = args,
                                method  = method)
            results.append(res)
        return results

#TODO: add multiprocessing
class MultiprocessingSolver(Solver):
    def run(self, args, method='MSE'):
        raise NotImplementedError


def solve(problemsGenerator, args, solver='sequential', method='MSE'):
    solver = {
        'sequential':   SequentialSolver,
        'MaxFlow':      MultiprocessingSolver
    }[solver](problemsGenerator)
    res = solver.run(args=args, method=method)
    return res
