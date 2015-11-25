library(MassTodon)
library(dplyr)
# ?MassTodonize

precursor 	<- MassTodon::ubiquitin
charge 		<- 10L
maxProtonsNo<- charge
cutIntensity<- 1500
cT 			<- 'C13H18N2O3'
# nT 			<- 'C13H18N2O3'
tolerance 	<- .05
path 		<- '~/Dropbox/Science/MassSpectrometry/MassTodonTests/workflow/2015_11/data/Melphalan_UBQ_100001000.txt'
cutIntensity<- 1500
cutPercentage<- 0
spectrum 	<- NA
ridge 		<- .1
multiplierQ <- 1.25
solver 		<- "kernlab"
parallel 	<- FALSE
cv 			<- FALSE
lambdas 	<- 10

Theory <- theory$new(
	precursorSequence 	= precursor,
	maxProtonsNo 		= charge,
	probabilityChebyshev= .01,
	cT = cT
)

Theory$formularize()

# Theory$formulas$forFit$monoisotopicMass %>% summary
# Theory$formulas$forFit[ which.max(Theory$formulas$forFit$monoisotopicMass), ]
# Theory <- theory$new(
# 	precursorSequence 	= precursor,
# 	maxProtonsNo 		= charge,
# 	probabilityChebyshev= .01,
# 	nT = nT
# )
# Theory$formularize()
# Theory$formulas$forFit$monoisotopicMass %>% summary
# Theory$formulas$forFit[ which.max(Theory$formulas$forFit$monoisotopicMass), ]

Theory$getMzStatistics()
Theory$getClusterEndings()

Empiria	<- empiriaStics$new(
	tolerance		= tolerance,
	path 			= path,
	Spectrum 		= spectrum,
	cutIntensity 	= cutIntensity,
	cutPercentage 	= cutPercentage
)

Empiria$preprocessSpectrum()
Empiria$trimEmpiricalData()
Fit <- fitSticks$new(
	formulas 		= Theory$formulas$forFit,
	presentElements = Theory$presentElements,
	clustersEndings = Theory$clustersEndings,
	spectrum  		= Empiria$spectrum$forFit,
	ridge 			= ridge,
	multiplierQ		= multiplierQ,
	parallel 		= parallel,
	tolerance 		= tolerance,
	solver 			= solver,
	lambdas 		= lambdas,
	cv      		= cv
)

Fit$findSpectrumOutsideClusters()
Fit$matchFormulasSpectrum()
Fit$prepareSubproblems()

####################################
# Adding the bloody new criterion for subproblem inclusion.
subproblem 			<- sapply(Fit$subproblems, function(x) x$empiria$intensity %>% sum) %>% which.max %>% Fit$subproblems[[.]]

	# Problems with 29 and 38
subproblem <- Fit$subproblems[[29]]

clusterIntensity 	<- sum( subproblem$empiria$intensity )
presentElements 	<- Fit$presentElements
theory 				<- subproblem$theory[, c("active",presentElements), drop = FALSE ]
theory$formulaName 	<- row.names( theory )
formula <- theory[1,]

# Error in apply(assignment[mostImportantPeaks, ], 1, any) : 
#   dim(X) must have a positive length

peaksForAssignment2 <- apply(
	theory, 1,
	function( formula ){
			# chooses the correct atoms' counts
		atomCount 	<- lapply(as.list( formula[presentElements] ), as.integer)
		aCharge 	<- as.integer(formula['active'])
		formulaName <- as.character(formula['formulaName'])
	
		brain <- BRAIN::useBRAIN( aC = atomCount, stopOption = 'coverage', nrPeaks = 10000, coverage = .9999)
	
		mzTheory<- (brain$masses)/aCharge
		
		theoreticalDistibution <- data.frame(
			mzTheory 	= as.numeric(mzTheory),
			probability	= brain$isoDistr,
			formulaNo  	= formulaName
		)

		distrWithOptimalSet <- theoreticalDistibution %>% dplyr::arrange( desc(probability) ) %>% mutate(
			totProb = cumsum(probability)
		) 

		i = 1
		stop = FALSE
		while( i <= nrow(distrWithOptimalSet) & !stop ){
			if( distrWithOptimalSet$totProb[i] > .8) stop = TRUE
			i = i+1
		}
		i <- min(i-1, nrow(distrWithOptimalSet))
		criticalSet <- distrWithOptimalSet %>% slice(1:i)

		if( nrow(criticalSet) < 3 ){
			criticalSet <- distrWithOptimalSet[1:min(3, nrow(distrWithOptimalSet)),]
		}

		mostImportantPeaks <- theoreticalDistibution$mzTheory %in%  criticalSet$mzTheory

		assignment 	<- apply(
			subproblem$empiria, 1,
			function( sample ) 	
				sample[1] + tolerance >= theoreticalDistibution$mzTheory &
				sample[1] - tolerance <= theoreticalDistibution$mzTheory
		)		 

		allCriticalPeaksPresent <- all( apply( assignment[mostImportantPeaks,,drop=FALSE], 1, any ) )		

		if( allCriticalPeaksPresent ){
			
			empiricalPeaksWithinRange <- apply( assignment, 2, any )

			if( sum( empiricalPeaksWithinRange ) > 1L ){

				maxEmpiricalCharge <- 
					( 1/min( 
						diff( subproblem$empiria$mz[empiricalPeaksWithinRange] 
					) ) ) %>%
					round %>% as.integer
					
		
				# if( maxEmpiricalCharge == formula[1] )	return(
				if( maxEmpiricalCharge >= formula['active'] ) return(
					list(
						theoreticalDistibution 	= theoreticalDistibution,
						assignment 				= assignment
					)
				) 
			} 

		}
		
		return( NA )
	}
)

nonZeroCoefficients2 <- !sapply(
	peaksForAssignment2,
	function( peaks ) class(peaks) == "logical"
)

nonZeroCoefficientsNo2 	<- sum( nonZeroCoefficients2 )
peaksForAssignment2	<- peaksForAssignment2[ nonZeroCoefficients2 ]

is.error <- function(x) inherits(x, "try-error")
res <- sapply(
		Fit$subproblems,
		function(subproblem){		
			clusterIntensity 	<- sum( subproblem$empiria$intensity )
			presentElements 	<- Fit$presentElements
			theory 				<- subproblem$theory[, c("active",presentElements), drop = FALSE ]
			theory$formulaName 	<- row.names( theory )
			
			peaksForAssignment2 <- try(
				apply(
				theory, 1,
				function( formula ){
						# chooses the correct atoms' counts
					atomCount 	<- lapply(as.list( formula[presentElements] ), as.integer)
					aCharge 	<- as.integer(formula['active'])
					formulaName <- as.character(formula['formulaName'])
					brain 	<- BRAIN::useBRAIN( aC = atomCount, stopOption = 'coverage', nrPeaks = 10000, coverage = .9999)
					mzTheory<- (brain$masses)/aCharge
					theoreticalDistibution <- data.frame(
						mzTheory 	= as.numeric(mzTheory),
						probability	= brain$isoDistr,
						formulaNo  	= formulaName
					)
					distrWithOptimalSet <- theoreticalDistibution %>% dplyr::arrange( desc(probability) ) %>% mutate(
						totProb = cumsum(probability)
					) 
					i = 1
					stop = FALSE
					while( i <= nrow(distrWithOptimalSet) & !stop ){
						if( distrWithOptimalSet$totProb[i] > .8) stop = TRUE
						i = i+1
					}
					i <- min(i-1, nrow(distrWithOptimalSet))
					criticalSet <- distrWithOptimalSet %>% slice(1:i)
					if( nrow(criticalSet) < 3 ){
						criticalSet <- distrWithOptimalSet[1:min(3, nrow(distrWithOptimalSet)),]
					}
					mostImportantPeaks <- theoreticalDistibution$mzTheory %in%  criticalSet$mzTheory
					assignment 	<- apply(
						subproblem$empiria, 1,
						function( sample ) 	
							sample[1] + tolerance >= theoreticalDistibution$mzTheory &
							sample[1] - tolerance <= theoreticalDistibution$mzTheory
					)		 
					allCriticalPeaksPresent <- all( apply( assignment[mostImportantPeaks,,drop=FALSE], 1, any ) )
					if( allCriticalPeaksPresent ){
						empiricalPeaksWithinRange <- apply( assignment, 2, any )
						if( sum( empiricalPeaksWithinRange ) > 1L ){
							maxEmpiricalCharge <- 
								( 1/min( 
									diff( subproblem$empiria$mz[empiricalPeaksWithinRange] 
								) ) ) %>%
								round %>% as.integer
							# if( maxEmpiricalCharge == formula[1] )	return(
							if( maxEmpiricalCharge >= formula['active'] ) return(
								list(
									theoreticalDistibution 	= theoreticalDistibution,
									assignment 				= assignment
								)
							) 
						} 
					}		
					return( NA )
				}
			)
		)
	
		errorIs <- is.error(peaksForAssignment2)
		return(errorIs)
	}
)

res[res]

library(ggplot2)
reduced2 <- bind_rows(
	peaksForAssignment2 %>%
	lapply(
		function(x) x$theoreticalDistibution %>% data.frame
	) %>% dplyr::bind_rows() %>%
	rename(mz = mzTheory) %>%
	mutate( tag= 'theory', intensity = -200000*probability ) %>% select(-probability),
	subproblem$empiria %>% mutate( tag = 'empiria' )
)


reduced <- bind_rows(
	peaksForAssignment %>%
	lapply(
		function(x) x$theoreticalDistibution %>% data.frame
	) %>% dplyr::bind_rows() %>%
	rename(mz = mzTheory) %>%
	mutate( tag= 'theory', intensity = -200000*probability ) %>% select(-probability),
	subproblem$empiria %>% mutate( tag = 'empiria' )
)

methodsComparison <- 
	bind_rows(
		reduced %>% mutate( method = 'original' ),
		reduced2%>% mutate( method = 'optimal80%' )
	) %>% ggplot(
		aes(
			x = mz,
			y = intensity,
			colour = tag	
		)
	) +
	geom_segment(
		aes(
			xend 	= mz,
			yend 	= 0
		)
	) + 
	facet_grid( method ~ .)

pdf( 
	file = '/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodonTests/devel/tooManyChargesFit/comparison.pdf',
	width 	= 100, 
	height 	= 10
)	
	print( methodsComparison )
dev.off()


subproblem$empiria %>% tbl_df %>% plot(type='h')














####################################

Fit$getFit()
Fit$getEstimates()
Fit$divideSpectrum()

clustersEndings <- Theory$clustersEndings
Spectrum <- c(
	Fit$spectrum,
	Empiria$spectrum
)

aggregateIntensities 	<- list(
	antePreprocessing 	= Empiria$totalIntensity,
	noSmallPeaks 		= Empiria$totalIntensityTrimmed
)

aggregateIntensities 	<- c(
	aggregateIntensities,
	lapply(
		Spectrum,
		function(x) { sum(x$intensity) }
	)
)

aggregateIntensities <-	unlist(aggregateIntensities)
aggregateIntensities/aggregateIntensities[2]