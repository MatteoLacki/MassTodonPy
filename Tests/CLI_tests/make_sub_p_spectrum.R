setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/CLI_tests')
spec = read.csv('subP_spectrum.csv')
plot(spec, type='h')
write.table(spec,'subP_spectrum.txt', row.names = FALSE, col.names = FALSE)
