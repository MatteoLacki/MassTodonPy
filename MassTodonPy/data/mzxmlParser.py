import 	rpy2.robjects as r

r.r.source('mzxmlParser.R', verbose=False)


try:
	mzXMLlibrary = importr('readMzXmlData')
except NameError:
	print('Stay awake - you have to choose the CRAN mirror.')
	utils = importr('utils')
	utils.install_packages('readMzXmlData')
	mzXMLlibrary = importr('readMzXmlData')

readMzXml  	= rO.r['readMzXmlFile']
unlist 		= rO.r['unlist']

def openMzXML(file, getMetaData=True ):
	spectraR = readMzXml(	mzXmlFile 		= file,
							removeMetaData 	= not getMetaData,
							verbose 		= False 			)
	spectra = []
	metaInfo= []
	for spectrum, metaData in spectraR:
		spectra.append([(mass, intensity) for mass, intensity in zip(spectrum[0],spectrum[1])])
		metaData = unlist(metaData)
		metaInfo.append( dict( zip(metaData.names, metaData) ) )
	return spectra, metaInfo
