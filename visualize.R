library(dplyr)
library(ggplot2)

D <- read.table('test.csv', header=TRUE)

sum(D[,1])
qplot(D[,1], bins=1000)
