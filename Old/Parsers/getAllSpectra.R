suppressPackageStartupMessages(require(readMzXmlData))

getRspectra <- function(path)
    readMzXmlFile(path, removeMetaData=TRUE, verbose = FALSE)
# x = getRspectra()
# str(x)
