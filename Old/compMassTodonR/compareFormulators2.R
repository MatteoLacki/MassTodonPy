library(MassTodonR)
library(dplyr)
library(rPython)

source('formularize.R')
# fasta   <- ubiquitin
fasta   <- 'RPKPQQFFGLM'
Q       <- 3L

Theory 	<- MassTodonize(fasta, Q)
formulas<- Theory$formulas$forFit %>% tbl_df

pasteFormula <- Vectorize(function(C,H,N,O,S){
    x <- c(C=C,H=H,N=N,O=O,S=S)
    x <- x[x>0]
    paste0(names(x),x,collapse="",sep="")
})

formulas_R = formulas %>% tbl_df %>%
            mutate(form=pasteFormula(C,H-active-neutral,N,O,S)) %>%
            rename(type=fragmentName, q=active, g=neutral) %>%
            select(type, form, q, g) %>%
            mutate(where='R', key=paste0(form,q,g))

formulas_R %>% filter(type=='c2')

formulas_py = read.csv("/Users/matteo/Documents/MassTodon/MassTodonPy/Old/compMassTodonR/ubiquitin_formulas_py.csv")
formulas_py = formulas_py %>% mutate(where='Py', key=paste0(form,q,g))


formulas_py %>% filter(!key %in% formulas_R$key )
formulas_R %>% filter(!key %in% formulas_py$key )

formulas_py %>% filter(type=='c5')
formulas_R  %>% filter(type=='c5')


# full_join(formulas_py, formulas_R, by=c(key='key')) %>% data.frame

# formulas_R %>% filter(type=='precursor')
# nrow(formulas_py)
# nrow(formulas_R)



# setdiff(unique(formulas_R$key), unique(formulas_py$key))
