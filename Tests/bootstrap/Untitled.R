library(jsonlite)
library(tidyverse)

x = fromJSON('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/cvxopt_optim/error_structures')

spectrum = tibble(
  mz = x[[2]],
  intensity = x[[3]]
)

hahaha = function(m){
  r = as.data.frame(m)
  colnames(r) = c('mz', 'intensity','tag')
  return(r)
}

D = lapply(x[[1]], hahaha) %>% bind_rows()
D = bind_rows(
  D, spectrum %>% mutate(tag = 'real')
)

D2 = D %>% group_by(tag) %>% mutate( prob = intensity/sum(intensity))


D2 %>% ggplot(aes(x=mz, y=prob, color=tag)) + geom_point()
