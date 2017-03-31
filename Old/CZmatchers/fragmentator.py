# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki.
#
#   This file is part of czMatchmaker.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   czMatchmaker is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with czMatchmaker.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

import numpy as np
import scipy.optimize as sciopt

def collect_fragments(fragments, Q):
    cNodes = []; zNodes = []
    for _,cz,q,I in fragments[['tag','active','estimates']].itertuples():
        if cz[0] == 'c':
            cNodes.append((cz,q,I))
        else:
            zNodes.append((cz,q,I))
    nodes = []
    intensities = []
    interactions= []
    costs = []
    for cz, q, I in cNodes + zNodes:
    intensities.append(I)
        nodes.append( (cz,q) )
        interactions.append(frozenset([(cz,q)]))
        costs.append( Q-1-q )
    for c, cQ, cI in cNodes:
        for z, zQ, zI in zNodes:
            if cQ + zQ < Q:
                interactions.append(frozenset([(c,cQ),(z,zQ)]))
                costs.append(Q-1-cQ-zQ)
    A = np.zeros((len(nodes), len(interactions)))
    for i,n in enumerate(nodes):
        for j, interaction in enumerate(interactions):
            if n in interaction:
                A[i,j] = 1

    print costs, A, intensities
    optim_result = sciopt.linprog( costs, A_eq=A, b_eq=intensities, options={"disp": False})
    results = {}
    for inter, I in zip(interactions, optim_result.x):
        if I > 0:
            results[inter] = I
    return results, optim_result
