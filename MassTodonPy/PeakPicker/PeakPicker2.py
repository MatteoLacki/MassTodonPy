#
# -*- coding: utf-8 -*-
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

from collections import Counter
from collections import defaultdict
from intervaltree import Interval as II
from intervaltree import IntervalTree as InTree
import networkx as nx
from six.moves import zip
from time import time

inf = float('inf')


class PeakPicker(object):
    """Class for peak picking."""

    def __init__(self,
                 spectrum,
                 molecules,
                 isotopic_calculator,
                 mz_precision_digits=0.05,
                 _verbose=True):
        """Initialize peak picker."""
        self.spectrum = spectrum
        self.molecules = molecules
        self.isotopic_calculator = isotopic_calculator
        self.mz_precision_digits = mz_precision_digits
        self.M = 0  # the number of molecule nodes
        self.I = 0  # the number of isotopologue nodes
        self._verbose = _verbose

    def represent_as_Graph(self):
        """Prepare the Graph based on mass spectrum and the formulas.

        Returns
        -------
        Graph : graph
            The deconvolution graph.

        """
        prec = self.mz_precision_digits
        Exps = InTree(II(mz - prec, mz + prec, (mz, intensity))
                      for mz, intensity in self.spectrum)
