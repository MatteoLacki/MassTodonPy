setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')


library('ggplot2')
library('tidyverse')

D = read.csv("shortData.csv") %>% tbl_df
R = read.csv("substanceP.csv") %>% tbl_df


p = D %>% select(mz_L, mz_R, tot_estimate, tot_intensity) %>%
  rename( L = mz_L, R=mz_R ) %>%
	mutate( mean_mz = factor((L+R)/2) ) %>%
	gather( "tag", "value", 3:4 ) %>%
	ggplot()+
	geom_bar(aes(x=mean_mz, y=value, fill=tag), stat='identity', position = 'dodge') +
	theme_dark()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	geom_hline(yintercept = 100)

# library('plotly')
# ggplotly(p)

ggsave(filename='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/spectra_plots/test_plot.pdf', plot=p, width=40, height=6)


# filter(mean_mz > 1340, mean_mz < 1360)  %>%
D2= D %>% mutate( mean_mz = (L+R)/2) %>% 	
	gather( "tag", "value", 1:2 ) %>%
	mutate( mean_mz = factor(mean_mz)) 

	
D2_I = D2 %>% filter(tag=='I')
D2_E = D2 %>% filter(tag=='E')
	

p =	ggplot() +
	geom_bar(data=D2_I, aes(x=mean_mz, y=value), stat = "identity", fill='white') + 
  	geom_bar(data=D2_E, aes(x=mean_mz, y=value), stat = "identity", fill='red') + 
	theme_dark()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	geom_hline(yintercept = 100)

ggplotly(p)

W =	D2 %>%
	ggplot( aes(x=mean_mz, y=value, color=tag)) +
	geom_point() +
	geom_linerange(aes(x=mean_mz, ymin = 0, ymax = value))

ggplotly(W)


q = R %>% mutate(mz = factor(mz)) %>%
	ggplot() + 
	geom_bar(aes(mz, intensity), stat='identity', col='green') +
	theme_dark()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

R %>% 
	ggplot() + 
	geom_point(aes(mz, intensity), col='green') +
	theme_dark()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


R2 = R %>% filter(mz > 1340, mz < 1360) 
plot(R2$mz, R2$intensity, type='h')
	




# p = D %>% ggplot() +
# 	geom_rect(aes(xmin=L, xmax=R, ymin=0, ymax=I)) +
# 	theme_minimal()
# # p = data_frame(x = 1:100, y = rnorm(100)) %>%
# # 	ggplot(aes(x=x,y=y))+
# # 	geom_point()

ggplotly(p)

