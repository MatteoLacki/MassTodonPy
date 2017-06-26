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
from    MassTodonPy.Misc import cvxopt_wrapper
from    time import time
from    multiprocessing import Pool
from    itertools import repeat
from    collections import Counter

# class Solver(object):
#     def __init__(self, problems, verbose=False):
#         self.problems = problems
#         self.verbose  = verbose
#         self.stats    = Counter()
#
#     def run(self, args, method='MSE', max_times_solve=5):
#         raise NotImplementedError
#
#
# class SequentialSolver(Solver):
#     def run(self,
#             args,
#             method='MSE',
#             max_times_solve=5
#         ):
#
#         results = []
#
#         with cvxopt_wrapper():
#             for SG in self.problems:
#                 i = 0
#                 stop = False
#                 T0 = time()
#
#                 while not stop:
#                     T00 = time()
#                     res = deconvolve(   SG      = SG,
#                                         args    = args,
#                                         method  = method)
#                     T01 = time()
#                     if self.verbose:
#                         print 'CVXOPT call lasted', T01-T00
#
#                     i += 1
#                     if res['status'] == 'optimal':
#                         stop = True
#                     if i == max_times_solve:
#                         stop = True
#                         print 'Deconvolution proved non optimal', max_times_solve, 'times'
#                 results.append(res)
#                 T1 = time()
#                 self.stats['Deconvolution Total T'] = T1-T0
#
#                 if self.verbose:
#                     print 'Solved problem in', T1-T0, 'It was big as', len(SG)
#                     print
#
#         return results

def worker(worker_args):
    SG, args, method, max_times_solve, verbose = worker_args
    i = 0
    stop = False
    T0 = time()
    while not stop:
        T00 = time()
        res = deconvolve(   SG      = SG,
                            args    = args,
                            method  = method )
        T01 = time()
        if verbose:
            print 'CVXOPT call lasted', T01-T00
        if res['status'] == 'optimal':
            stop = True
        if i == max_times_solve:
            stop = True
            print 'Deconvolution proved non optimal', max_times_solve, 'times'
    T1 = time()
    if verbose:
        print 'Solved problem in' , T1-T0, 'It was big as', len(SG)
    return res

#
# class MultiprocessingSolver(Solver):
#     def run(self,
#             args,
#             method,
#             max_times_solve = 5,
#             multiprocesses_No = None
#         ):
#         '''Run the multiprocesses solver.'''
#             #Start solving bigger, i.e. graphs with more nodes, problems first
#         self.problems.sort(reverse=True, key=len) # len = len(SG) = #Nodes
#         pool_args = zip(self.problems,
#                         repeat(args),
#                         repeat(method),
#                         repeat(max_times_solve),
#                         repeat(self.verbose) )
#
#         with cvxopt_wrapper():
#             T0 = time()
#             P = Pool(multiprocesses_No)
#             results = P.map(    worker,
#                                 pool_args    )
#             P.close()
#             P.join()
#             T1 = time()
#
#         self.stats['Deconvolution Total T'] = T1-T0
#
#         if self.verbose:
#             print 'Solved problem in', T1-T0
#             print
#
#         return results

def solve(  problems,
            args,
            solver              = 'sequential',
            multiprocesses_No   = None,
            method              = 'MSE',
            max_times_solve     = 5,
            verbose             = False    ):

    '''Wrapper over the solver class.

    Runs the solver with a given set of inputs.'''

    stats = Counter()

    with cvxopt_wrapper():
        T0 = time()
        if solver == 'sequential':
            results = [worker((SG, args, method, max_times_solve, verbose)) for SG in problems]
        elif solver == 'multiprocessing':
            #Start solving bigger, i.e. graphs with more nodes, problems first
            problems.sort(reverse=True, key=len) # len = len(SG) = #Nodes
            pool_args = zip(self.problems,
                            repeat(args),
                            repeat(method),
                            repeat(max_times_solve),
                            repeat(self.verbose) )
            P = Pool(multiprocesses_No)
            results = P.map( worker, pool_args )
            P.close()
            P.join()
        else:
            print "'solver' should have been either 'sequential' or 'multiprocessing'. Stick to the name conventions, please."
            raise NotImplementedError
        T1 = time()

    stats['Deconvolution Total T'] = T1-T0
    if verbose:
        print 'Solved problem in', T1-T0
        print

    return res
