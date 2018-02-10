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
from collections import Counter
from cvxopt import matrix, spmatrix, sparse, spdiag, solvers, setseed
from future.builtins import super
from math import sqrt
import networkx as nx
from random import randint

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Deconvolution.Misc import diag, normalize_rows

#TODO: try to eliminate copying while instantiating
#TODO: turn the matrix generators into iterators
#      so to use the with different solvers easily
class DeconvolutionProblem(nx.Graph):
    """Prepare and solve one deconvolution problem."""
    def __init__(self,
                 data=None,
                 L1_flow=0.001,
                 L2_flow=0.001,
                 L1_intensity=0.001,
                 L2_intensity=0.001,
                 max_times=10,
                 show_progress=False,
                 maxiters=1000,
                 **kwds):
        """Deconvolve the isotopic signal using the tolerance-interval-method.

        Parameters
        ==========
        data : nx.Graph inputs
            Necessary for compatibility with 'networkx.Graph'.
        L1_flow : float
            L1 penalty for high flows of intensities.
        L2_flow : float
            L2 penalty (a.k.a. ridge regression like) for high flows of intensities.
        L1_intensity : float
            L1 penalty for high intensity estimates.
        L2_intensities : float
            L2 penalty (a.k.a. ridge regression like) for high intensities.
        max_times : int
            The maximal number of times to run CVXOPT.
        show_progress : boolean
            Show progress of the CVXOPT calculations.
        maxiters : int
            Maximum number of iterations for the CVXOPT algorithm.
        kwds
            Arguments for other functions.
        """
        global solvers
        super().__init__(data, **kwds)
        self.count_nodes_and_edges()
        self.max_times = int(max_times)
        self.L1_flow = float(L1_flow)
        self.L2_flow = float(L2_flow)
        self.L1_intensity = float(L1_intensity)
        self.L2_intensity = float(L2_intensity)
        self.get_P_q()
        self.get_A_b()
        self.get_G_h()
        self.get_initvals()
        solvers.options['show_progress'] = bool(show_progress)
        solvers.options['maxiters'] = int(maxiters)
        self.solve()
        self._get_global_fit_quality()

    def count_nodes_and_edges(self):
        """Tag nodes and edges with distinct id numbers."""
        cnts = Counter()
        for N in self:
            self.node[N]['cnt'] = cnts[N[0]]
            cnts[N[0]] += 1
            if N[0] is not 'I':  # N-I in { M-I, G-I }
                for I in self[N]:
                    self[N][I]['cnt'] = cnts[N[0] + I[0]]
                    cnts[N[0] + I[0]] += 1
        cnts['var'] = cnts['GI'] + cnts['M']
        self.__dict__.update({k + str('_no'): v
                              for k, v in cnts.items()})

    def node_iter(self, node_type, *args):
        """Iterate over the nodes of a given type."""
        if args:
            for N, N_data in self.nodes(data=True):
                if N[0] is node_type:
                    out = [N]
                    out.extend(N_data.get(a, None) for a in args)
                    yield tuple(out)
        else:
            for N in self:
                if N[0] is node_type:
                    yield N

    def get_P_q(self):
        """Get the cost function parameters: matrix P and vector q.

        Notes
        =====
        0.5 <x|P|x> + <q|x> + L1_flow * sum x + L2_flow * sum x^2 + L1_intensity * sum alpha + L2_intensity * sum alpha^2
        """
        q_list = []
        P_list = []
        for G, intensity in self.node_iter('G', 'intensity'):
            degree = self.degree(G)
                                            # L1 penalty for x
            q_list.append(matrix(-intensity + self.L1_flow,
                                 size=(degree, 1)))
            ones = matrix(1.0, (degree, 1))
                                      # L2 penalty for x
            P_g  = ones * ones.T + diag(self.L2_flow, degree)
            P_list.append(P_g)
                           # L1 penalty for alphas
        q_list.append(matrix(self.L1_intensity, (self.M_no, 1)))
        self.q = matrix(q_list)
                         # L2 penalty for alphas
        P_list.append(diag(self.L2_intensity, self.M_no))
        self.P = spdiag(P_list)

    def get_initvals(self):
        """Initial values of flows."""
        self.initvals= {}
        self.initvals['x'] = matrix(0.0, (self.var_no, 1))

    def get_G_h(self):
        """Prepare for conditions Gx <= h"""
        self.G = diag(-1.0, self.var_no)
        self.h = matrix(0.0, size=(self.var_no, 1))

    def get_A_b(self):
        """Prepare for conditions Ax = b."""
        A_x=[]; A_i=[]; A_j=[]
        for M, M_cnt in self.node_iter('M','cnt'):
            for I in self[M]:
                I_cnt = self.node[I]['cnt']
                A_x.append(-self.node[I]['probability'])
                A_i.append(I_cnt)
                A_j.append(M_cnt + self.GI_no)
                for G in self[I]:
                    if not G == M:
                        A_x.append(1.0)
                        A_i.append(I_cnt)
                        A_j.append(self[G][I]['cnt'])
        self.A = spmatrix(A_x, A_i, A_j, size=(self.I_no, self.var_no))
        self.A = normalize_rows(self.A)
        self.b = matrix(0.0, (self.I_no, 1))

    def sum_intensities(self):
        """Sum real intensities in a graph."""
        return sum(i for G, i in self.node_iter('G', 'intensity'))

    def solve(self):
        """Solve the given deconvolution problem."""
        stop = False
        solved = False
        iteration = 1
        print(self.max_times)
        while not stop or iteration <= self.max_times:
            setseed(randint(0, 1000000))
            try:
                self.solution = solvers.qp(self.P, self.q, self.G,
                                           self.h, self.A, self.b,
                                           initvals=self.initvals)
                if self.solution['status'] is 'optimal':
                    stop = True
                    solved = True
            except ValueError as e:
                solved = False
                print('iteration')
            iteration += 1
        if not solved:
            print("Tried {} times and no optimum reached.".format(self.max_times))

        # filling up the graph
        X = self.solution['x']
        self.min_mz = infinity
        for M, M_cnt in self.node_iter('M', 'cnt'):
            self.node[M]['estimate'] = X[self.GI_no + M_cnt]
            self.node[M]['molecule'].intensity = X[self.GI_no + M_cnt]
        for G in self.node_iter('G'):
            self.node[G]['estimate'] = 0.0
            try:
                min_mz = self.node[G]['min_mz']
            except KeyError:
                min_mz = self.node[G]['mz']
            self.min_mz = min(self.min_mz, min_mz)
            for I in self[G]:
                idx = self[G][I]['cnt']
                self.node[G]['estimate'] += X[idx]
                self[G][I]['estimate'] = X[idx]


    def _get_global_fit_quality(self):
        self.status = self.solution['status']
        if self.status is 'optimal':
            self.L1_error = sum(abs(i - e) for G, i, e in
                                self.node_iter('G',
                                               'intensity',
                                               'estimate'))
            self.L2_error = sqrt(sum((i - e)**2 for G, i, e in
                                     self.node_iter('G',
                                                    'intensity',
                                                    'estimate')))
            self.underestimates = sum(max(e - i, 0) for G, i, e in
                                      self.node_iter('G',
                                                     'intensity',
                                                     'estimate'))
            self.overestimates = sum(max(i - e, 0) for G, i, e in
                                     self.node_iter('G',
                                                    'intensity',
                                                    'estimate'))
        else:
            self.overestimates = 0.0
            self.L1_error = self.sum_intensities()
            self.L2_error = sqrt(sum(i**2 for G, i in
                                     self.node_iter('G', 'intensity')))
            self.underestimates = self.L1_error

        self.total_intensity = sum(i for _, i in self.node_iter('G', 'intensity'))

    def global_fit_quality(self):
        return {'l1': self.L1_error,
                'l2': self.L2_error,
                'overestimates': self.overestimates,
                'underestimates': self.underestimates,
                'status': self.status,
                'total intensity': self.total_intensity}
