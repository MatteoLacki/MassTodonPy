library(IsoSpecR)
suppressMessages(library(dplyr))
library(ggplot2)

mol <- c(C=100, H=200, N=20)
prob<- .999
Q   <- 6

data(isotopicData)
H1 <- isotopicData$IsoSpec %>% filter(isotope == 'H1')
H1mass <- H1$mass
result <- IsoSpecify(mol, prob, fancy=TRUE)
result <- result %>% mutate( mass = mass + Q*H1mass, H1 = H1 + Q, type='proxy')

mol2 <- mol
mol2['H']<- mol2['H']+Q
result2 <- IsoSpecify(mol2, prob, fancy=TRUE) %>% mutate(type='real')

# bind_rows(  result %>% select(mass, prob) %>% mutate(type='proxy'),
#             result2%>% select(mass, prob) %>% mutate(type='real')) %>%
#     ggplot(aes(x=factor(mass),y=prob)) +
#     geom_bar(stat='identity', aes(fill = type), position = "dodge")+
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#     xlab('mass')

bind_rows(  result %>% select(mass, prob),
            result2%>% select(mass, prob) ) %>%
    ggplot(aes(x=mass,y=prob, color=type, shape=type)) +
    geom_point()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab('mass')


IsoSpecify( c(H=20), prob, fancy=TRUE) # This is the maximal

result
