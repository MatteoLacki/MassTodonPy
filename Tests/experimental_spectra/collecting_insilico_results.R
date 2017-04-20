suppressPackageStartupMessages(library("jsonlite"))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

path='/Users/matteo/Documents/MassTodon/in_silico_results/wloczykij_latest.json'
R <- read_json(path, simpifyVector=TRUE)
P <- function(x, l=1) str(x, max.level=l )
P(R)
# R <- tbl_df(read.csv('/Users/matteo/Documents/MassTodon/in_silico_results/analysis.csv'))

Probs <-
	c(  0.99999943,  0.99919849,  0.98831003,  0.95667342,  0.90814916,
        0.8518448 ,  0.79452072,  0.73981509,  0.68929589,  0.64340622,
        0.60204112,  0.56485823,  0.53143491,  0.50134371,  0.4741858 ,
        0.44960344,  0.4272823 ,  0.40694922,  0.388368  ,  0.37133488,
        0.35567401,  0.34123345,  0.32788163,  0.31550433,  0.30400205,
        0.29328792,  0.28328577,  0.27392868,  0.2651576 ,  0.25692031,
        0.24917047,  0.24186683,  0.2349726 ,  0.22845482,  0.22228395,
        0.21643338,  0.21087914,  0.20559956,  0.20057498,  0.19578757,
        0.19122111,  0.18686079,  0.18269311,  0.17870568,  0.17488715,
        0.17122712,  0.16771598,  0.16434491,  0.16110575,  0.15799095,
        0.15499353,  0.15210704,  0.14932545,  0.14664319,  0.14405505,
        0.14155621,  0.13914213,  0.1368086 ,  0.13455168,  0.13236767 )

Sigma2Prob <- data.frame(sigma=sort(unique(R$sigma)), prob=Probs)
R <- left_join(R, Sigma2Prob, by='sigma')

labelNice <-
	function(molsNo) ifelse( molsNo==1,
		paste0(molsNo, ' MOLECULE'),
		paste0(molsNo, ' MOLECULES'))

PlotWrap <-
	R %>%
	filter(molsNo<7) %>%
	ggplot(aes(x=1-prob,y=error, color=x0))+
	geom_point()+
	theme_minimal()+
	scale_x_continuous(limits = c(0, 1), labels = scales::percent)+
	scale_y_continuous(limits = c(0, 1))+
	scale_colour_gradient(trans = "log10")+
	scale_color_continuous(name="Ions Number")+
	geom_smooth(method='lm',formula=y~x,color='red', se=FALSE)+
	facet_wrap(~labelNice(molsNo))+
	ylab('Error')+
	xlab('Probability of a Molecule not being within a Tolerance Interval')

PlotWrap
ggsave('/Users/matteo/Documents/MassTodon/in_silico_results/PlotWrap.pdf', plot=PlotWrap )


PlotOne <-
	R %>%
	filter(molsNo<7) %>%
	ggplot(aes(x=1-prob,y=error))+
	geom_point(size=.6, aes(color=factor(molsNo)), alpha=.5)+
	theme_minimal()+
	scale_x_continuous(limits = c(0, 0.88), labels = scales::percent)+
	scale_y_continuous(limits = c(0, 1))+
	# scale_alpha_gradient(trans = "log10")+
	# scale_alpha_continuous(name="Ions Number")+
	geom_smooth(method='lm', se=FALSE, aes(group=molsNo,color=factor(molsNo)))+
	ylab('Error')+
	coord_fixed()+
	xlab('Probability of a Molecule not being within a Tolerance Interval')


PlotOne
ggsave('/Users/matteo/Documents/MassTodon/in_silico_results/PlotOne.pdf', plot=PlotOne )
