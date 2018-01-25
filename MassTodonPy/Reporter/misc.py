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

from bokeh.palettes import viridis, Colorblind, Paired, Set3, Set1
from itertools import cycle


class Cluster(object):
    __slots__ = ('mz_L', 'mz_R', 'mz_left', 'mz_right', 'intensity',
                 'intensity_d', 'estimate', 'estimate_d')

    def __init__(self, **kwds):
        for s in self.__slots__:
            setattr(self, s, kwds.get(s, None))

    def __repr__(self):
        res = []
        for s in self.__slots__:
            v = getattr(self, s)
            res.append( "{0}={1}".format(s,v))
        return "Cluster(" + ", ".join(res) + ")"


class PeakGroup(object):
    """A class representing the outcomes for one experimental grouping."""
    __slots__ = ('mz_L', 'mz_R', 'mz_left', 'mz_right', 'intensity',
                 'intensity_d', 'estimate', 'estimate_d')

    def __init__(self, **kwds):
        for s in self.__slots__:
            setattr(self, s, kwds.get(s, None))

    def __repr__(self):
        res = []
        for s in self.__slots__:
            v = getattr(self, s)
            res.append( "{0}={1}".format(s,v))
        return "PeakGroup(" + ", ".join(res) + ")"


class Brick(object):
    __slots__ = ('peak_group', 'top', 'bottom', 'color', 'intensity', 'molecule')

    def __init__(self, **kwds):
        for s in self.__slots__:
            setattr(self, s, kwds.get(s, None))

    def __repr__(self):
        res = []
        for s in self.__slots__:
            v = getattr(self, s)
            res.append( "{0}={1}".format(s,v))
        return "Brick(" + ", ".join(res) + ")"


class ColorGenerator(object):
    """Generator of colors with memoization."""

    def __init__(self, pallette='Set1'):
        self.color = cycle({'Colorblind': Colorblind[8],
                            'Paired': Paired[12],
                            'Set1': Set1[9],
                            'Set3': Set3[12]}[pallette])
        self.colors = {}

    def __call__(self, molecule):
        """Get the next color in sequence."""
        try:
            return self.colors[molecule]
        except KeyError:
            self.colors[molecule] = next(self.color)
            return self.colors[molecule]


def make_string_represenation(name, digits):
    return "@" + str(name) + "{0." + "".join(["0"] * digits) + "}"
