from aminoAcid2 import AminoAcid

aas = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

for aa in aas:
	AminoAcid(aa).plot(target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/AAs/'+aa+'.pdf')

try:
	for aa in aas:
		X = AminoAcid(aa)
		X.add_OH()
		X.add_H()
		X.plot(target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/AAs/'+aa+'_with_OH_and_H.pdf')
except ValueError:
	
AminoAcid('W').plot()

x = AminoAcid('W')
x.plot()
x.add_OH()
x.add_H()
x.plot()


y = AminoAcid('P')
y.plot()

if x.G['name'] == 'Proline':
	x.G.add_vertex(name='HN1', elem='H')
	x.G.add_edge('N1','HN1')
else:
	x.G.vs.find('HNalpha')['name'] = 'HNalpha1'
	x.G.add_vertex(name='HNalpha2', elem='H')
	x.G.add_edge('Nalpha','HNalpha2') # Might this be another by bond? Assume not