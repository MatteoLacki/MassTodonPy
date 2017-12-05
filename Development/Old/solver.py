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
from MassTodonPy.Misc import cvxopt_wrapper
from time import time
from multiprocessing import Pool
from itertools import repeat
from collections import Counter


class Error_in_update_scaling(Exception):
    pass

def worker(worker_args):
    SG, args, method, max_times_solve, verbose = worker_args
    solved = False
    try:
        i = 0
        stop = False
        while not stop:
            i += 1
            res = deconvolve(   SG      = SG,
                                args    = args,
                                method  = method )
            if res['status'] == 'optimal':
                stop    = True
                solved  = True
            else:
                # for penalty in ('L1_x', 'L2_x', 'L1_alpha', 'L2_alpha'):
                #     args[penalty] *= 10.0
                if verbose:
                    print args
            if i == max_times_solve:
                stop = True
                print 'Deconvolution proved non optimal', max_times_solve, 'times'
    except Error_in_update_scaling as e:
        pass
    if verbose:
        if solved:
            print 'Solved problem. It was big as', len(SG)
        else:
            print 'Problem not solved. It was big as', len(SG)
    return res


def solve(  problems,
            args,
            solver              ='sequential',
            multiprocesses_No   = None,
            method              ='MSE',
            max_times_solve     = 5,
            verbose             = False):
    stats = Counter()

    with cvxopt_wrapper():
        T0 = time()
        if solver == 'sequential':
            results = [worker((SG, args, method, max_times_solve, verbose))
                       for SG in problems]
        elif solver == 'multiprocessing':
            # Start solving bigger, i.e. graphs with more nodes, problems first
            problems.sort(reverse=True, key=len) # len(SG) = #Nodes
            pool_args = zip(problems,
                            repeat(args),
                            repeat(method),
                            repeat(max_times_solve),
                            repeat(verbose) )
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

    stats.update(Counter([ r['status'] for r in results ]))

    return results, stats
