from formulaGenerator.fragments import makeFragments, ubiquitin
	

	# Making pure fragments:
makeFragments('AAAPPP')
makeFragments('A', fragmentTypes=['ax','cz'])
makeFragments('AAAA', fragmentTypes=['ax','cz'], innerFragments=True)
res = makeFragments(ubiquitin, innerFragments=True)
len(res)