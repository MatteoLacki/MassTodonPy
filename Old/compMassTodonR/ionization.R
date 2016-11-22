library(MassTodonR)
library(dplyr)
source('formularize.R')


fasta 	<- ubiquitin
ions 	<- formularize(fasta, proteinogenicAAs = MassTodonR::proteinogenicAAs)
maxProtonsNo <- 9L


####################### protonations
precursorProtonation<- do.call(
	rbind,
	lapply(
		1:maxProtonsNo,
		function( activeProtonation )
			data.frame(
				active 	= activeProtonation,
				neutral	= 0:(maxProtonsNo - activeProtonation)
			)
	)
)

formulas <- list()
formulas$forFit <- merge(
	do.call(
		rbind,
		apply(
			ions[ c('fragmentName', 'subsequence')],
			1,
			function(ion)
				data.frame(
					fragmentName = ion[1],
					switch( seqinr::s2c(ion[1])[1],
						'p' = precursorProtonation,
						subset(
							precursorProtonation,
							active <= (nchar(ion[2]) %/% 5) + 1 # The charge restriction
						)
					),
					row.names = NULL
				)
		)
	),
	ions
)

ionsType <- substring( formulas$forFit$fragmentName, 1, 1 )
cond <- ionsType == 'p' | (ionsType != 'p' & formulas$forFit$active < maxProtonsNo)

formulas$forFit <- formulas$forFit[cond,]




	### Adding the monoisotopic masses.
presentMonoIsotopes <- monoisotopes[ monoisotopes$element %in% presentElements, ]
rownames(presentMonoIsotopes) <- presentMonoIsotopes$element
		# reordering according to the order of atoms in presentElements
presentMonoIsotopes <- presentMonoIsotopes[presentElements,]
rownames(presentMonoIsotopes) <- NULL

formulas$forFit <<- formulas$forFit %>% 
	dplyr::mutate( H = H + active + neutral ) 

formulas$forFit <<- formulas$forFit %>% dplyr::mutate(
	monoisotopicMass = as.vector(as.matrix(formulas$forFit[,presentElements]) %*% presentMonoIsotopes$mass)
)	
