library(dplyr)
sh = function(x, ml=1) str(x, max.level = ml)

sh(data4plots)
sh(data4plots$WH)
sh(data4plots$WH$MassTodons$'FRL-010513-SUBP-WH000-WV300'$cTerminus)


sapply(data4plots$WH$MassTodons, function(x) x$nTerminus )
sapply(data4plots$WH$MassTodons, function(x) x$cTerminus )
sapply(data4plots$WV$MassTodons, function(x) x$cTerminus )
sapply(data4plots$WV$MassTodons, function(x) x$nTerminus )
