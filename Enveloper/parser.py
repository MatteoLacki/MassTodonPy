ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

def parse( precursor ):
	revPrecursor 	= precursor[::-1]
	fragments 		= [revPrecursor[0]]
	prolineChunk 	= revPrecursor[0] == 'P'
	for i in range( 1, len(revPrecursor)): #LGFFQQPKPR
		AA = revPrecursor[i]
		if  AA == 'P':		
			if prolineChunk:
				fragments[-1] += AA 			
			else:
				fragments.append('P')

			prolineChunk = True	
		else: 
			if prolineChunk:
				fragments[-1] += AA
			else:
				fragments.append(AA)		
			prolineChunk = False		
	if precursor[0] != 'P':
		fragments.append('')
	fragments = [ x[::-1] for x in fragments[::-1] ]
	residues  = [ len(x) for x in fragments ]
	return ( fragments, residues )	

# fragments, residues = parse(ubiquitin)
# print fragments