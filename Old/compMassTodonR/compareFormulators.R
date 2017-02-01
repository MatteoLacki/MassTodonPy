library(MassTodonR)
library(dplyr)
library(rPython)

source('formularize.R')

fasta   <- ubiquitin
Theory 	<- theory$new( precursorSequence = fasta, maxProtonsNo= 9L, probabilityChebyshev = .01 )
Theory$formularize()
formulas<- Theory$formulas$forFit

formulas <- formulas %>% tbl_df %>%
            mutate(form=paste0('C',C,'H',H,'N',N,'O',O,'S',S)) %>%
            rename(type=fragmentName, q=active, g=neutral) %>%
            select(type, form, q, g)


read.csv("/Users/matteo/Documents/MassTodon/MassTodonPy/Old/compMassTodonR/ubiquitin_formulas_py.csv")






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
# python.load("/Users/matteo/Documents/Science/MassTodon/MassTodonPy/formulaGenerator/callFromR.py")
#
#
# ### Checking if all fragments have the same formulas. # Perfect match.
# fragments <- python.get("fragments")
# fragments <- lapply(fragments, function(x) data.frame(name=x[[5]], t(x[[4]])) )
# fragments <- do.call(rbind,fragments) %>% tbl_df
# fragments <- fragments %>% arrange( name )
# fragments <- fragments %>% mutate( name=as.character(name) ) %>% arrange( name )
#
#
#
# fragments 	%>% filter(name == 'precursor')
# ions 		%>% filter(name == 'precursor')
#
# ions 		<- ions %>% arrange( name )
# fragments 	<- fragments %>% arrange( name)
# identical(ions, fragments %>% filter(name != 'c0')) # Perfect.
