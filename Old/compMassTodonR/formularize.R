add_vectors <- function(...) {
  a <- list(...)
  nms <- sort(unique(unlist(lapply(a, names))))
  out <- numeric(length(nms))
  names(out) <- nms
  for (v in a) out[names(v)] <- out[names(v)] + v

  out
}

getSubsequence 		<- function(
	fragmentName,
	precursor
){
	if( fragmentName == 'precursor') return( seqinr::c2s(precursor) )
	else {
		fragmentNameLength 	<- nchar( fragmentName )
		type 				<- substring( fragmentName, 1L, 1L)
		number  		 	<- as.integer(
			substring( fragmentName, 2L, fragmentNameLength)
		)

		sequenceLength 	<- length(precursor)

		switch(
			type,
			'c'= {
				if( number == 0L ) return('') else return(
					paste( precursor[ 1:number ], collapse="")
				)
			},
			'z'= return(
				paste( precursor[ (sequenceLength-number+1):sequenceLength ], collapse="" )
			),
			error("Not good naming convention while providing fragment name.
				Should be 'xN', where x is either 'c' or 'z' and N is a number.")
		)
	}
}

load('/Users/matteo/Documents/MassTodon/MassTodonR/R/sysdata.rda')
# fasta = MassTodonR::ubiquitin
formularize = function(){


	aminoAcidNo <- length( precursor )
	notProline 	<- !precursor == "P"

		# Preparing the c-z ions reparametrisation.
			# allAminoAcids are stated in counts of b-y ions.

		# Restrict to elements that are forming the molecule.
	aminoAcidResiduals <- lapply(precursor, function(aa) proteinogenicAAs[[aa]] )

		# Adding the modifications
	if( 'aaMods' %in% names(modifications) ){

		AAwithDiffIdx <- modifications$aaMods %>% dplyr::select( aminoAcidNo ) %>% unlist

		aminoAcidResiduals[ AAwithDiffIdx ] <-
			Map(
				function( AAatomCnt, diff ) add_vectors( unlist(diff), unlist(AAatomCnt) ),
				aminoAcidResiduals[AAwithDiffIdx],
				modifications$aaMods %>% dplyr::select(-aminoAcidNo ) %>% split(1:nrow(modifications$aaMods))
			)
	}

		# adding the c0 part.
			# Simply: add H for the non-residual n-Terminus H
			# Add the NH as part of the c0 fragment (reparametrization).
	c0 <- c( H = 2L, N = 1L )

	if( 'cTermMod' %in% names(modifications) ){
		c0 <- add_vectors( c0, unlist(modifications$cTermMod) )
	}

	z1idx <- length(aminoAcidResiduals)
	z1 <- aminoAcidResiduals[[ z1idx ]]

	z1 <- z1 %>%
		add_vectors( c( H =-1L,  N =-1L ) ) %>% # The NH cannot be stolen from the next AA, cause there is no next AA.
		add_vectors( c( H = 1L,  O = 1L ) )		# The residual form lacks OH.

	if( 'nTermMod' %in% names(modifications) ) z1 <-
		z1 %>% add_vectors( unlist(modifications$nTermMod)	)

	aminoAcidResiduals[[z1idx]] <- z1
	aminoAcidResiduals <- c(c0=list(c0), aminoAcidResiduals)

	AAR4precursorCalc <- bind_rows(aminoAcidResiduals %>% lapply(function(x) x %>% t %>% data.frame)) %>% as.matrix
	AAR4precursorCalc[ is.na(AAR4precursorCalc) ] <- 0L

			# Making precursors
	precursors <- cbind(
		as.data.frame(t(colSums(AAR4precursorCalc))),
		fragmentName 	= 'precursor',
		subsequence 	= seqinr::c2s( precursor )
	)
	rownames(precursors) <- 'precursor'

		# Here we prepare a modified list of residuals. We add an extra H
		# to the c0 fragment, assuring that the extra H gets counted in every
		# c fragment, and omitted in every z fragment.

	AAR4fragmentCalc <- AAR4precursorCalc
	if( !('H' %in% colnames(AAR4fragmentCalc)) ) stop("You're protein has no hydrogens left. If this is not an issue find this place in code and comment it.")

	AAR4fragmentCalc[1,'H'] <- AAR4fragmentCalc[1,'H'] + 1L

	if( any(AAR4fragmentCalc < 0) ) stop('Your modifications resulted in negative counts of atoms (told ya to be careful).')

		# These are all the elements that make up the
	presentElements <<- colnames(AAR4fragmentCalc)

		# This internal function will calculate the fragments
		# atomic compositions.
	getFragments <- function(
		ionType
	){
		fragments <- apply(
			AAR4fragmentCalc[
				if( ionType == 'c' )
					-nrow(AAR4fragmentCalc)	# we neglect z1
				else
					nrow(AAR4fragmentCalc):2 # we neglect c0
				,
			],
			2,
			cumsum
		)

		rownames(fragments) <- fragmentNames <- paste0(
			ionType,
			1:nrow(fragments) +
				if( ionType == 'c' ) -1L else 0L
		)

		fragments <- cbind(
			as.data.frame(fragments),
			fragmentName = fragmentNames
		)

		return(fragments)
	}


		# Making fragments
		#
		# observation by Frederik: indeed, we change the definition of the
		# c ions by substracting one hydrogen, for it seems more intuitive.
	cFragments <- getFragments('c')
	zFragments <- getFragments('z')


	# If there are prolines certain fragments won't appear.
	if( any(!notProline) )
	{
		prolineNoInSequence <- (1:aminoAcidNo)[!notProline]

			# No c0 or any c_{K-1} when P is the K-th aminoacid.
		cFragmentsAbsent <- paste0( 'c', c( 0, prolineNoInSequence - 1))

			# No z_{N-K+1} when P is the the K-th aminoacid in a N aminoacid long sequence.
		zFragmentsAbsent <- paste0( 'z', aminoAcidNo - prolineNoInSequence + 1)

		cFragments <- cFragments[ setdiff( rownames(cFragments), cFragmentsAbsent), ]
		zFragments <- zFragments[ setdiff( rownames(zFragments), zFragmentsAbsent), ]
	}


		# Adding subsequence tags.
	cFragments$subsequence <- sapply(
		as.character(cFragments$fragmentName),
		getSubsequence
	)

	zFragments$subsequence <- sapply(
		as.character(zFragments$fragmentName),
		getSubsequence
	)

	ions <- rbind(
		precursors,
		cFragments,
		zFragments
	)

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

	# formulas$forFit <<- merge(
	# 	do.call(
	# 		rbind,
	# 		apply(
	# 			ions[ c('fragmentName', 'subsequence')],
	# 			1,
	# 			function(ion)
	# 				data.frame(
	# 					fragmentName = ion[1],
	# 					switch( seqinr::s2c(ion[1])[1],
	# 						'p' = precursorProtonation,
	# 						'z' = subset(
	# 							precursorProtonation,
	# 							active <= (nchar(ion[2]) %/% 5) + 1 # The charge restriction
	# 						),
	# 						'c' = subset(
	# 							precursorProtonation %>% mutate( neutral = neutral-1L ),
	# 							active <= (nchar(ion[2]) %/% 5) + 1 # The charge restriction
	# 						)
	# 					),
	# 					row.names = NULL
	# 				)
	# 		)
	# 	),
	# 	ions
	# )

	potentialCharges <- function(ions, distanceBetweenCharges=5){
		aasNo 		= nchar(ion[2])
		potentialQ 	= aasNo %/% distanceBetweenCharges	
		if( aasNo %% distanceBetweenCharges > 0 ) potentialQ = potentialQ + 1
		potentialQ
	}

	formulas$forFit <<- merge(
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
							'z' = subset(
								precursorProtonation,
								active <= potentialCharges(ion) # The charge restriction
							),
							'c' = subset(
								precursorProtonation %>% mutate( neutral = neutral-1L ),
								active <= potentialCharges(ion) # The charge restriction
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

	formulas$forFit <<- formulas$forFit[cond,]

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

	return(invisible())
}
