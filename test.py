from formulaGenerator.fragments import makeFragments, ubiquitin, roepstorffy, fasta2atomCount, protonatedFragments

	# Getting atomic composition of a fasta.
ubi_atomCnt = fasta2atomCount(ubiquitin)
print(ubi_atomCnt)		


	# Making pure fragments:
# res = makeFragments('AAAPPP')
# res = makeFragments('A', fragmentTypes=['ax','cz'])
# res = makeFragments('AAAA', fragmentTypes=['ax','cz'], innerFragments=True)
res = makeFragments(ubiquitin, fragmentTypes=['cz'])
len(res) # Don't know what to kill? Kill runtime.


	
	# Fragments with protonation numbers.
res = list( protonatedFragments( ubiquitin, 20) )