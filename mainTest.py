%load_ext autoreload
%autoreload
from Formulator.makeSuperAtoms import cz_fragments, pandizeSubstances
from Formulator.protonations import protonate


ubiquitin   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
substanceP  = 'RPKPQQFFGLM'
fastas      = [substanceP, ubiquitin, ubiquitin+ubiquitin+ubiquitin]

Q = 10
list(protonate(Q,'c'))


# modifications = {   ('N',2) :       lCnt({'H': -1, 'O': +2, 'N': +3}),
#                     ('Calpha',2) :  lCnt({'H': -1, 'O': +2, 'N': +3}),
#                     ('Calpha',5) :  lCnt({'H': -2, 'S': +2, 'N': +2}),
#                     ('C',6) :       lCnt({'H': -2, 'S': +2, 'N': +2}) }

# modifications = {   ('N',2) :       lCnt({'H': -1, 'O': +2, 'N': +3}),
                    # ('Calpha',2) :  lCnt({'H': -1, 'O': +2, 'N': +3}),
                    # ('Calpha',5) :  lCnt({'H': -2, 'S': +2, 'N': +2}),
                    # ('C',6) :       lCnt({'H': -2, 'S': +2, 'N': +2}) }
modifications = {}
precursor, cFrags, zFrags = cz_fragments(substanceP, modifications)
pandizeSubstances(precursor, cFrags, zFrags)



# def sideChainsNo(fragment):
#     '''Finds the number of side chains on a fragment.'''
#     pos, leftAA, rightAA, _ = fragment
#     leftPos, rightPos       = pos
#     res = rightAA - leftAA + 1
#     if leftPos == 'R':
#         res -= 1
#     if rightPos== 'L':
#         res -= 1
#     if res < 0:
#         res = 0
#     return res


# def protonatedFragments(    fasta,
#                             max_charge,
#                             amino_acids_per_charge  = 5,
#                             fragmentTypes           = ['cz'],
#                             innerFragments          = False     ):
#     fragments = makeFragments(fasta, fragmentTypes, innerFragments)
#     for fragment in fragments:
#         maxQonFrag = round(sideChainsNo(fragment)/float(amino_acids_per_charge))
#         for q,g in getProtonation( max_charge, maxQonFrag ):
#             yield (q,g, fragment)

# res = list( protonatedFragments( ubiquitin, 20) )
# len(res)
