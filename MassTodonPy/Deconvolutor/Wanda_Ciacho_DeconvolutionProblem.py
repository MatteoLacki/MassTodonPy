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

from MassTodonPy.Deconvolutor.DeconvolutionProblem import DeconvolutionProblem


def l2_dot_product(mean_0, mean_1, var_0, var_1):
    """Get the L2 dot product of two Gaussians."""
    out = log(pi) + log(var_0 + var_1) - (mean_0 - mean_1)**2/(var_0 + var_1)
    return exp(-.5*out)

class DeconvolutionProblem(DeconvolutionProblem):
    def __init__(self,
                 data=None,  
                 L1=0.001,
                 L2=0.001,
                 max_times=10,
                 show_progress=False,
                 maxiters=1000,
                 **attr):
        super(DeconvolutionProblem, self).__init__(data, **attr)
        self.L1 = L1
        self.L2 = L2
        self.max_times = max_times
        solvers.options['show_progress'] = show_progress
        solvers.options['maxiters'] = maxiters
        self.count_nodes_and_edges()
        self.get_P_q()
        self.get_A_b()
        self.get_G_h()
        self.get_initvals()
        
    def get_P_q(self):
        """Get the cost function parameters: matrix P and vector q.

        Notes
        =====
        0.5 <x|P|x> + <q|x> + L1_x * sum x + L2_x * sum x^2 + L1_alpha * sum alpha + L2_alpha * sum alpha^2
        """
        pass

    def get_A_b(self):
        pass

    def get_G_h(self):
        pass
