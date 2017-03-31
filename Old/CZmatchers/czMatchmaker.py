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

import pandas as pd
from dplython import DplyFrame, X, sift, select, mutate, rename, summarize, group_by, DelayFunction
from fragmentator import collect_fragments
from collections import Counter
from misc import crossprod, Round

def czMatchmaker(data , Q, precursor_fasta):
    data = pd.read_csv("/Users/matteo/Documents/czMatchmaker/data/examplaryData.csv")
    data = DplyFrame(data)
    precursors = data >> \
    	sift( X.tag == 'precursor' ) >> \
    	select( X.active, X.neutral, X.estimates)

    fragments = data >> sift( X.tag != 'precursor' ) >> \
    	group_by( X.tag, X.active, X.broken_bond ) >> \
    	summarize( estimates = X.estimates.sum() )

    I_on_fragments 	= {}
    optiminfos 		= {}
    for break_point, data in fragments.groupby('broken_bond'):
        pairing, optiminfo 			= collect_fragments(data, Q)
        I_on_fragments[break_point] = pairing
        optiminfos[break_point] 	= optiminfo

    cations_fragmented_I = sum(
        sum( I_on_fragments[bP][p] for p in I_on_fragments[bP])
        for bP in I_on_fragments )

    I_no_reactions = precursors >> \
        sift( X.active==Q, X.neutral == 0) >> \
        select( X.estimates )

    I_no_reactions = I_no_reactions.values.flatten()[0]

    prec_ETnoD_PTR_I = precursors >> \
        sift( X.active != Q ) >> \
        rename( ETnoD  = X.neutral, I = X.estimates ) >> \
        mutate( PTR    = Q - X.ETnoD - X.active ) >> \
        select( X.ETnoD, X.PTR, X.I )

    I_prec_no_frag = prec_ETnoD_PTR_I >> \
        summarize( I = X.I.sum() )

    I_prec_no_frag = I_prec_no_frag.values.flatten()[0]

    precursorNoReactions = precursors >> \
        sift( X.active == Q ) >> \
        select( X.estimates )

    prec_ETnoD_PTR_I = prec_ETnoD_PTR_I >> mutate(
            I_PTR 	= crossprod(X.PTR, X.I), \
            I_ETnoD = crossprod(X.ETnoD, X.I) ) >> \
        summarize( I_PTR = X.I_PTR.sum(), I_ETnoD = X.I_ETnoD.sum() )

    I_PTR_no_frag, I_ETnoD_no_frag = prec_ETnoD_PTR_I.values.flatten()

    prob_PTR   = I_PTR_no_frag/(I_PTR_no_frag + I_ETnoD_no_frag)
    prob_ETnoD = 1. - prob_PTR

    I_frags = dict(
        ( bP, sum( I_on_fragments[bP][pairing] for pairing in I_on_fragments[bP] ) )
        for bP in I_on_fragments 	)

    I_frag_total = sum( I_frags[bP] for bP in I_frags )

    prob_frag = Counter(dict( (int(bP), I_frags[bP]/I_frag_total) for bP in I_frags))
    prob_frag = [ prob_frag[i] for i in range(len(precursor_fasta)) ]

    I_frags_PTRETnoD_total = sum(
        ( Q-1-sum( q for cz, q in pairing ) ) * I_on_fragments[bP][pairing] for bP in I_on_fragments for pairing in I_on_fragments[bP] 	)

    anion_meets_cation = I_frags_PTRETnoD_total + I_PTR_no_frag + I_ETnoD_no_frag
    prob_fragmentation = I_frags_PTRETnoD_total / anion_meets_cation
    prob_no_fragmentation = 1 - prob_fragmentation

    prob_no_reaction = I_no_reactions/( I_no_reactions + I_frag_total + I_prec_no_frag )
    prob_reaction = 1. - prob_no_reaction

    res = {}
    res['reaction'] = (prob_reaction, prob_no_reaction)
    res['fragmentation'] = (prob_fragmentation, prob_no_fragmentation)
    res['fragmentation_amino_acids'] = tuple(prob_frag)
    return res
