


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


# i = IsoSpecPy.IsoSpec.IsoFromFormula("H2O1", 0.9)
# print i.getConfs()


####################################
