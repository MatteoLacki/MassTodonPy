library(MassTodonR)
library(dplyr)
library(rPython)

source('formularize.R')

fasta 		<- ubiquitin
# fasta 		<- 'AAA'

proteinogenicAAs <- MassTodonR::proteinogenicAAs
proteinogenicAAs[['O']] <- NULL
proteinogenicAAs[['U']] <- NULL
finalFastaTest <- names(proteinogenicAAs) %>% seqinr::c2s()
finalFastaTest <- c(finalFastaTest,finalFastaTest,finalFastaTest) %>% seqinr::c2s()
ions <- formularize(finalFastaTest)
ions <- ions %>% tbl_df %>% select(-subsequence) %>% rename( name=fragmentName ) %>% transmute(
		name=name, H=H, C=C, S=S, O=O, N=N ) %>% arrange( name )
ions <- ions %>% mutate( name=as.character(name) ) %>% arrange( name )





### Comparing Residuals # Perfect match
# python.load("/Users/matteo/Documents/Science/MassTodon/MassTodonPy/formulaGenerator/callFromR.py")
# proteinogenicAAs
# AAs <- python.get("aminoAcids")
# AAs <- AAs %>% lapply(function(x) as.data.frame(t(x))) %>% bind_rows(.id="aminoAcids")
# AAs[is.na(AAs)] <- 0L

# AAsR<- proteinogenicAAs %>% lapply(function(x) as.data.frame(t(x))) %>% bind_rows(.id="aminoAcids")
# AAsR[is.na(AAsR)] <- 0L

# AAs <- AAs %>% transmute(aminoAcids, C,H,N,O,S)
# AAsR <- AAsR %>% transmute(aminoAcids, C,H,N,O,S,Se)

# AAs
# AAsR

# full_join(AAs, AAsR, by='aminoAcids') # Perfect match


# Load/run the main Python script
python.load("/Users/matteo/Documents/Science/MassTodon/MassTodonPy/formulaGenerator/callFromR.py")


### Checking if all fragments have the same formulas. # Perfect match.
fragments <- python.get("fragments")
fragments <- lapply(fragments, function(x) data.frame(name=x[[5]], t(x[[4]])) )
fragments <- do.call(rbind,fragments) %>% tbl_df
fragments <- fragments %>% arrange( name )
fragments <- fragments %>% mutate( name=as.character(name) ) %>% arrange( name )



fragments 	%>% filter(name == 'precursor')
ions 		%>% filter(name == 'precursor')

ions 		<- ions %>% arrange( name )
fragments 	<- fragments %>% arrange( name) 
identical(ions, fragments %>% filter(name != 'c0')) # Perfect.


