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


def parse_single_modification(mod):
    modification = {}
    for elem_no in mod.split(','):
        elem, no = elem_no.split(':')
        no = int(no)
        modification[elem] = no
    return modification


def parse_plain_config_file(path):
    '''Parse a text file in a format similar to:

    fasta = RPKPQQFFGLM
    precursor_charge = 3
    modification C11 = H:1, O:-1, N:1
    cut_off = 100
    mz_prec = .05
    '''
    config_tmp = {}
    with open(path, 'r') as f:
        for line in f:
            key, value = re.sub('[\s+]|[\n]','',line).split('=')
            config_tmp[key] = value
    config = {}
    for k in config_tmp:
        if 'modification' == k[0:12]:
            value = parse_single_modification(config_tmp[k])
            k = k[12:]
        else:
            try:
                value = float(config_tmp[k])
            except ValueError:
                value = config_tmp[k]
        config[k] = value
    return config

parse_plain_config_file('example_config_plain_text.txt')
