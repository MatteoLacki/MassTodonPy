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

def get_subsequence(fasta, name, modifications):
    '''For purpose of ETDetective, to bridge gaps in definitions.'''
    if modifications == {'C11': {'H': 1, 'N': 1, 'O': -1}}:
        suffix = '*'
    else:
        suffix = ''
    if name[0]=='p':
        return '*'+fasta+suffix
    if name[0]=='z':
        return fasta[ len(fasta)-int(name[1:]): ] + suffix
    else:
        return '*' + fasta[ 0:int(name[1:]) ]


def results_to_etdetective( raw_estimates, fasta, modifications={}, threshold = 0.0 ):
    '''Reorganize results of MassTodonPy to fit ETDetective format.

    Parameters
    ----------
    raw_estimates : list
        A list of MassTodon raw deconvolution outputs
    fasta: str
        The amino acid sequence.
    modifications: dict
        A dictionary containing the differences in atomic compositions between the target modified amino acid and the input amino acid.
    threshold : float
        Show outputs with estimates higher than that number.

    Returns
    -------
        A dictionary with keys corresponding to tuples (<sequence>, <charge>, <quenched charge>, <name of molecule>) and values corresponding to the estimate.
    '''

    etdetective_input = {}
    for subproblem in raw_estimates:
        for r in subproblem['alphas']:
            if r['estimate'] > threshold:
                seq = get_subsequence( fasta, r['molType'], modifications )
                etdetective_input[(seq, r['q'], r['g'], r['molType'])] =  r['estimate']
    return etdetective_input
