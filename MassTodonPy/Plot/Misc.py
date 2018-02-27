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

from collections import defaultdict

def aggregate_fragments(report, fasta, offset):
    """Aggregate the fragments for the plot."""
    data = defaultdict(list)
    for d in report.iter_aggregated_fragments_intensity():
        for k, v in d.items():
            data[k].append(v)
    data['x'] = list(range(len(fasta)))
    data['x_left'] = [i-offset for i in range(len(fasta))]
    data['x_right'] = [i+offset for i in range(len(fasta))]
    return data
