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

from MassTodonPy.Deconvolutor.Misc import diag, normalize_rows

class DeconvolutionProblem(nx.Graph):
    """Prepare and solve one deconvolution problem."""
    def __init__(self, 
                 data=None,
                 L1_x=0.001, 
                 L2_x=0.001,
                 L1_alpha=0.001,
                 L2_alpha=0.001,
                 max_times=10,
                 show_progress=False,
                 maxiters=1000,
                 **attr):
        global solvers
        super().__init__(data, **attr)
        self.count_nodes_and_edges()
        self.max_times = max_times
        self.L1_x = L1_x
        self.L2_x = L2_x
        self.L1_alpha = L1_alpha
        self.L2_alpha = L2_alpha
        self.get_P_q()
        self.get_A_b()
        self.get_G_h()
        self.get_initvals()
        solvers.options['show_progress'] = show_progress
        solvers.options['maxiters'] = maxiters

    def count_nodes_and_edges(self):
        """Add numbers to graph nodes and edges."""
        cnts = Counter()
        for N in self:
            self.node[N]['cnt'] = cnts[N[0]]
            cnts[N[0]] += 1
            if N[0] is not 'I':  # number edges
                for I in self[N]:
                    self[N][I]['cnt'] = cnts[N[0]+I[0]]
                    cnts[N[0]+I[0]] += 1
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
        0.5 <x|P|x> + <q|x> + L1_x * sum x + L2_x * sum x^2 + L1_alpha * sum alpha + L2_alpha * sum alpha^2
        """
        q_list = []
        P_list = []
        for G, intensity in self.node_iter('G', 'intensity'):
            degree = self.degree(G)
                                            # L1 penalty for x
            q_list.append(matrix(-intensity + self.L1_x,
                                 size=(degree, 1)))
            ones = matrix(1.0, (degree, 1))
                                      # L2 penalty for x
            P_g  = ones * ones.T + diag(self.L2_x, degree)
            P_list.append(P_g)
                           # L1 penalty for alphas
        q_list.append(matrix(self.L1_alpha, (self.M_no, 1)))
        self.q = matrix(q_list)
                         # L2 penalty for alphas
        P_list.append(diag(self.L2_alpha, self.M_no))
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

    def get_L1_error(self):
        """Get the L1 error.

        The sum of absolute values of differences between estimates and real intensities."""
        return sum(abs(i - e) for G, i, e in
                   self.node_iter('G', 'intensity', 'estimate'))

    def get_L2_error(self):
        """Get the L2 error.

        Get the Euclidean distance between estimates and real intensities."""
        return sqrt(sum((i - e)**2 for G, i, e in
                    self.node_iter('G', 'intensity', 'estimate')))

    def get_L1_signed_error(self, sign=1.0):
        """Get the L1 signed error.

        If sign=1, returns how much estimates are above real values.
        If sign=-1, returns how much estimates are below real values."""
        return sum(max(sign*(e - i), 0) for G, i, e in
                   self.node_iter('G', 'intensity', 'estimate'))

    def sum_intensities(self):
        """Sum real intensities in a graph."""
        return sum(i for G, i in self.node_iter('G', 'intensity'))

    def solve(self):
        """Solve the given deconvolution problem."""
        stop = False
        solved = False
        iteration = 1
        while not stop or iteration is self.max_times:
            setseed(randint(0, 1000000))
            self.solution = solvers.qp(self.P, self.q, self.G,
                                       self.h, self.A, self.b,
                                       initvals=self.initvals)
            if self.solution['status'] is 'optimal':
                stop = True
                solved = True
            iteration += 1
        if not solved:
            print("Tried {} times and no optimum reached.".format(self.max_times))

        # filling up the graph
        X = self.solution['x']
        self.alphas = []
        for M, M_cnt in self.node_iter('M', 'cnt'):
            self.node[M]['estimate'] = X[self.GI_no + M_cnt]
            self.alphas.append(self.node[M])
        for G in self.node_iter('G'):
            self.node[G]['estimate'] = 0.0
            for I in self[G]:
                idx = self[G][I]['cnt']
                self.node[G]['estimate'] += X[idx]
                self[G][I]['estimate'] = X[idx]

    def report(self):
        if self.solution['status'] is 'optimal':
            res = {'alphas': self.alphas,
                   'L1_error': self.get_L1_error(),
                   'L2_error': self.get_L2_error(),
                   'underestimates': self.get_L1_signed_error(sign=1.0),
                   'overestimates': self.get_L1_signed_error(sign=-1.0),
                   'status': self.solution['status']}
        else:
            res = {'status': 'ValueError',
                   'L1_error': self.sum_intensities(),
                   'overestimates': 0.0}
        return res
