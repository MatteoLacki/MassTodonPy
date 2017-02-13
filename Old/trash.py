
# pos = nx.shell_layout(P)
# # pos = nx.random_layout(P)
# nx.draw(P, pos=pos)
# plt.show()


# cc.node['I62585']     #  'I62585'
# cc['I62585']          # adjacency list of generalized isotopologue 'I62585'
# cc.nodes(data=True)   # all nodes with data
# Counter(contains_experimental_peaks(cc) for cc in ccs)
# good_ccs = [ cc for cc in ccs if contains_experimental_peaks(cc) ]
# bad_ccs  = [ cc for cc in ccs if not contains_experimental_peaks(cc) ]
# Counter(len(cc) for cc in good_ccs)
# Counter(len(cc.edges()) for cc in good_ccs)

# cc = good_ccs[12]
# len(cc)
# [ (len(cc), len(cc.edges())) for cc in good_ccs]
