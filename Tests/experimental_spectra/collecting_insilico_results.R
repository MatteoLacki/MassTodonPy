suppressPackageStartupMessages(library("jsonlite"))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(hash))

# path='/Users/matteo/Documents/MassTodon/in_silico_results/wloczykij_latest.json'
path='/Users/matteo/Documents/MassTodon/in_silico_results/uniform.json'
R <- read_json(path, simpifyVector=TRUE)
P <- function(x, l=1) str(x, max.level=l )

# sapply(R, function(x) c( sum(unlist(x[[2]][[1]])),sum(unlist(x[[2]][[2]])) ))
# R[[1]]

labelNice <-
	function(molsNo) ifelse( molsNo==1,
		paste0(molsNo, ' MOLECULE'),
		paste0(molsNo, ' MOLECULES'))


getDF = function(R, error_function){
    res = sapply( R,
        function(x) return(c(unlist(x[[1]]), error_function(unlist(x[[2]][[1]]), unlist(x[[2]][[2]]))))
    ) %>% t %>% tbl_df %>% rename(sigma=V1, prob=V2, molsNo=V3, x0=V4, error=V5) %>% mutate(molsNo=factor(molsNo))
    levels(res$molsNo) = labelNice(levels(res$molsNo))
    res
}

makePlotWrap = function(Res) Res %>%
	ggplot(aes(x=sigma,y=error))+
	geom_point(aes(color=x0))+
	theme_minimal()+
	scale_x_continuous(labels = scales::percent)+
	# scale_y_continuous(limits = c(0, .6), labels = scales::percent)+
	scale_colour_gradient(trans = "log10", name="Ions Number")+
	facet_wrap(~molsNo)+
	ylab('Error')+
	xlab('Probability of a Molecule not being within a Tolerance Interval')


ER0 = function(p,q) sum(abs(p-q))/2
makePlotWrap(getDF(R,ER0))

ER1 = function(p,q){
    lp = log(p)
    lq = log(q)
    sum(p*lp + q*lq - p*lq - q*lp)/(2*log(length(p)))
}
makePlotWrap(getDF(R,ER1))

ER2 = function(p,q) crossprod(p,q)
makePlotWrap(getDF(R,ER2))


ER3 = function(p,q) exp(log(crossprod(p,q)) - 0.5 * ( log(crossprod(p)) + log(crossprod(q)) ) )

ER3 = function(p,q) sqrt(1-(exp(log(crossprod(p,q)) - 0.5 * ( log(crossprod(p)) + log(crossprod(q)) ))^2 ))
plot_correl = makePlotWrap(getDF(R,ER3))

plot_correl + scale_y_continuous(labels = scales::percent) +
    geom_smooth(se=FALSE, color='red')+
    ylab('Correlation')


ER2 = function(p,q) sqrt(crossprod(p-q))
makePlotWrap(getDF(R,ER2))

pnorm = function(d)
    function(p,q) exp(log(sum(abs(p-q)**d))/d)
makePlotWrap(getDF(R,pnorm(1.5)))

ER3 = function(p,q) sum( abs(p-q)/(p+q) )/length(p)
makePlotWrap(getDF(R,ER3))

ER4 = function(p,q) max( abs(p-q)/(p+q) )
makePlotWrap(getDF(R,ER4))

ER45 = function(p,q) mean( abs(p-q)/(p+q) )
makePlotWrap(getDF(R,ER45))

ER5 = function(p,q) sum( abs(p-q)/p )/length(p)
makePlotWrap(getDF(R,ER5))

getDF(R,cor) %>% filter(error<0)




sapply( R,
    function(x) all(unlist(x[[2]][[1]]) >= 0.0001)
) %>% t %>% all








ggsave('/Users/matteo/Documents/MassTodon/in_silico_results/PlotWrap.pdf', plot=PlotWrap )


PlotOne <-
	Res0 %>%
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
