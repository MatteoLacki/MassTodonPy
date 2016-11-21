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

fasta = MassTodonR::ubiquitin
proteinogenicAAs = MassTodonR::proteinogenicAAs
formularize <- function(fasta, proteinogenicAAs = MassTodonR::proteinogenicAAs)
{
	precursor 	<- seqinr::s2c(fasta)
	aminoAcidNo <- length( precursor )
	notProline 	<- !precursor == "P"

		# Preparing the c-z ions reparametrisation.
			# allAminoAcids are stated in counts of b-y ions.

		# Restrict to elements that are forming the molecule.
	aminoAcidResiduals <- lapply(precursor, function(aa) proteinogenicAAs[[aa]] )

		# adding the c0 part.
			# Simply: add H for the non-residual n-Terminus H
			# Add the NH as part of the c0 fragment (reparametrization).
	c0 <- c( H = 2L, N = 1L )
	z1idx <- length(aminoAcidResiduals)
	z1 <- aminoAcidResiduals[[ z1idx ]] 
	z1 <- z1 %>%
		add_vectors( c( H =-1L,  N =-1L ) ) %>% # The NH cannot be stolen from the next AA, cause there is no next AA.
		add_vectors( c( H = 1L,  O = 1L ) )		# The residual form lacks OH. 

	aminoAcidResiduals[[z1idx]] <- z1
	aminoAcidResiduals <- c(c0=list(c0), aminoAcidResiduals)
	
	AAR4precursorCalc <- dplyr::bind_rows(aminoAcidResiduals %>% lapply(function(x) x %>% t %>% data.frame)) %>% as.matrix	
	AAR4precursorCalc[ is.na(AAR4precursorCalc) ] <- 0L 

			# Making precursors
	precursors <- cbind(
		as.data.frame(t(colSums(AAR4precursorCalc))),
		fragmentName 	= 'precursor',
		subsequence 	= seqinr::c2s( precursor )
	)
	rownames(precursors) <- 'precursor'
	
		# Here we prepare a modified list of residuals. We add an extra H
		# to the c1 fragment, assuring that the extra H gets counted in every
		# c fragment, and omitted in every z fragment. The extra H comes from the C=O turning into C-O-H with 
		# H coming from the ionized residual to the right. c0 hasn't got a free C=O to do this form of things.

	AAR4fragmentCalc <- AAR4precursorCalc
	if( !('H' %in% colnames(AAR4fragmentCalc)) ) stop("You're protein has no hydrogens left. If this is not an issue find this place in code and comment it.")
	
	AAR4fragmentCalc[1,'H'] <- AAR4fragmentCalc[1,'H'] + 1L # This line added an extra H.

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
		getSubsequence,
		precursor=precursor
	)

	zFragments$subsequence <- sapply(
		as.character(zFragments$fragmentName),
		getSubsequence,
		precursor=precursor
	)

	ions <- rbind(
		precursors,
		cFragments,
		zFragments
	)

	return(ions)
}