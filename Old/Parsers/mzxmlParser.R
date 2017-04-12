suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(readMzXmlData))

trim_spectrum <- function(spectrum, cut_off_intensity=50, round_mass_digits=2){
	spectrum 	<- spectrum$spectrum
	intensity 	<- spectrum$intensity
	mz 			<- spectrum$mass[ intensity > cut_off_intensity ]
	intensity 	<- intensity[ intensity > cut_off_intensity ]

	data_frame(mz=round(mz,round_mass_digits), intensity=intensity) %>%
	group_by(mz) %>%
	summarize(intensity=sum(intensity))
}

aggregate_spectrum <- function(
    path,
    cut_off_intensity = 50,
    round_mass_digits = 2
){
	D <- readMzXmlFile(path)
    lapply(D,
        trim_spectrum,
        cut_off_intensity,
        round_mass_digits
    ) %>% bind_rows %>%
	group_by(mz) %>%
    summarize(intensity=sum(intensity)) %>%
    as.list
}

# res <- aggregate_spectrum('data.mzXML')
# file = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/Ubiquitin_ETD_10 ms_1071.mzXML'
# aggregate_spectrum(file) %>% plot(type='h')
# D <- readMzXmlFile(file)
# length(D[[1]]$spectrum$intensity)
# lapply(D, trim_spectrum, cut_off_intensity=1000 )
