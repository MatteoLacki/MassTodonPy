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
from    MassTodonPy.Deconvolutor import deconvolve
from    time  import time
from    multiprocessing import Pool
from    itertools import repeat
from    collections import Counter

class Solver(object):
    def __init__(self, problemsGenerator, verbose=False):
        self.prob_gen = problemsGenerator
        self.verbose  = verbose
        self.stats    = Counter()

    def run(self, args, method='MSE', max_times_solve=5):
        raise NotImplementedError


class SequentialSolver(Solver):
    def run(self, args, method='MSE', max_times_solve=5):
        results = []
        for SG in self.prob_gen:
            i = 0
            stop = False
            T0 = time()

            while not stop:
                T00 = time()
                res = deconvolve(   SG      = SG,
                                    args    = args,
                                    method  = method)
                T01 = time()
                if self.verbose:
                    print 'CVXOPT call lasted', T01-T00

                i += 1
                if res['status'] == 'optimal':
                    stop = True
                if i == max_times_solve:
                    stop = True
                    print 'Deconvolution proved non optimal', max_times_solve, 'times'
            results.append(res)
            T1 = time()
            self.stats['Deconvolution Total T'] = T1-T0

            if self.verbose:
                print 'Solved problem in', T1-T0, 'It was big as', len(SG)
                print


        return results


def helper(helper_args):
    SG, args, method, verbose = helper_args
    T0 = time()
    res = deconvolve(   SG      = SG,
                        args    = args,
                        method  = method )
    if res['status'] != 'optimal':
        print 'Deconvolution proved non optimal',
    T1 = time()
    if verbose:
        print 'Solved problem in' , T1-T0, 'It was big as', len(SG)

    return res


class MultiprocessingSolver(Solver):
    def run(self, args, method, multiprocesses_No = None):
        T0 = time()

            #TODO make this multiprocessed too.
        problems = list(self.prob_gen)

            #Start solving bigger, i.e. graphs with more nodes, problems first
        problems.sort(reverse=True, key=len) # len = len(SG) = #Nodes

        pool_args = zip(    problems,
                            repeat(args),
                            repeat(method),
                            repeat(self.verbose) )

        P = Pool(multiprocesses_No)

        results = P.map(    helper,
                            pool_args    )

        P.close()
        P.join()
        T1 = time()
        self.stats['Deconvolution Total T'] = T1-T0

        if self.verbose:
            print 'Solved problem in', T1-T0
            print

        return results


def solve(  problemsGenerator,
            args,
            solver              = 'sequential',
            multiprocesses_No   = None,
            method              = 'MSE',
            max_times_solve     = 5,
            verbose             = False    ):
    '''Wrapper over the solver class.

    Runs the solver with a given set of inputs.'''

    if solver == 'sequential':
        S   = SequentialSolver(problemsGenerator, verbose)
        res = S.run(args   = args,
                    method = method,
                    max_times_solve = max_times_solve   )

    if solver == 'multiprocessing':
        S = MultiprocessingSolver(problemsGenerator, verbose)
        res = S.run(args   = args,
                    method = method,
                    multiprocesses_No = multiprocesses_No   )

    return res
