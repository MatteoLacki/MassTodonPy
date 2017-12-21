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
from cvxopt import matrix, spmatrix, sparse, spdiag, solvers, setseed
from math import sqrt, exp, log, pi
import networkx as nx

from MassTodonPy.Deconvolution.DeconvolutionProblem import DeconvolutionProblem


def l2_dot_prod(mean_0, mean_1, var_0, var_1):
    """Get the L2 dot product of two Gaussians."""
    out = log(2) + log(pi) + log(var_0 + var_1) - (mean_0 - mean_1)**2/(var_0 + var_1)
    return exp(-out/2.0)


class GaussianDeconvolutionProblem(DeconvolutionProblem):
    def __init__(self,
                 data=None,
                 L1=0.001,
                 L2=0.001,
                 max_times=10,
                 show_progress=False,
                 maxiters=1000,
                 sigma2=.1, 
                 ni2=.1,
                 **attr):
        # This is important to initiate a graph!!!!
        #   calls the grandparent class
        super(DeconvolutionProblem, self).__init__(data, **attr)
        self.L1 = L1
        self.L2 = L2
        self.sigma2 = sigma2
        self.ni2 = ni2
        self.max_times = max_times
        solvers.options['show_progress'] = show_progress
        solvers.options['maxiters'] = maxiters
        self.count_nodes_and_edges()
        self.get_P_q()
        self.get_G_h()
        self.get_initvals()

    def count_nodes_and_edges(self):
            """Tag nodes and edges with distinct id numbers."""
            self.M_no = 0
            for M in self:
                if M[0] is 'M':
                    self.node[M]['cnt'] = self.M_no
                    self.M_no += 1
            self.var_no = self.M_no


    def get_P_q(self):
        """Get the cost function parameters: matrix P and vector q.

        Notes
        =====
        0.5 <x|P|x> + <q|x> + L1_x * sum x + L2_x * sum x^2 + L1_alpha * sum alpha + L2_alpha * sum alpha^2
        """
        q = [0.0] * self.M_no
        for M in self.node_iter('M'):
            M_idx = self.node[M]['cnt']
            for I in self[M]:
                prob = self.node[I]['probability']
                I_mz = self.node[I]['mz']
                for E in self[I]:
                    if E[0] is 'E':  # automatically neglect M-I-∅
                        E_mz = self.node[E]['mz']
                        intensity = self.node[E]['intensity']
                        product = l2_dot_prod(I_mz,
                                              E_mz,
                                              self.sigma2,
                                              self.ni2)
                        q[M_idx] += prob * intensity * product
        self.q = matrix(q)

    def get_G_h(self):
        pass

    def solve(self):
        pass