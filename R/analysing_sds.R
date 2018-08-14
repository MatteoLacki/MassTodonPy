library(readr)
library(tidyverse)
library(MASS)
library(ggthemes)
D = read_csv("data/data_4_sd.csv")
D = D %>% mutate(log_mean = log(mean),  log_sd = log(sd), log_intensity=log(intensity)) 

M = lm(sd ~ mean + cnt + intensity, D)
summary(M)

# check the distribution of peak counts in different m/z
with(D, plot(mean, cnt))
D %>% ggplot(aes(x=ordered(cnt), y=mean)) + geom_boxplot() + theme_tufte() + coord_flip()
# ergo: no clear dependcence upon the m/z position and the number of counts.
table(D$cnt)
qplot(x=cnt, data=D, geom='bar') # 5-9 counts dominate.

# check all other things :D
plot(D, pch='.')
with(D, plot(mean, spread, col=cnt))
with(D, plot(sd, spread, col=cnt))
with(D, plot(mean, spread, cex=intensity))
with(D, plot(mean, spread))
qplot(x = mean, y = spread, size = sd, color = intensity, data = D)
with(D, plot(cnt, spread))
with(D %>% filter(cnt %in% 5:8), plot(sd, log_intensity, col=cnt))

D %>% 
  filter(cnt %in% 5:12,
         log_sd > -6) %>%
  qplot(x=spread, y=log_intensity, data=.) + 
  geom_smooth(method='rlm') +
  geom_smooth(method='lm', color='red') +
  facet_grid(.~cnt)



D %>% ggplot(aes(mean, spread, color=sd)) + geom_point()
D %>% 
  filter(cnt >= 5) %>%
  ggplot(aes(sd, spread, color=ordered(cnt))) +
  geom_point(size=.2) +
  geom_smooth(se = F, method='lm') + 
  theme_tufte()

table(D$cnt)

D %>% 
  filter(cnt == 5) %>%
  ggplot(aes(sd, sqrt(intensity))) +
  geom_point() +
  geom_smooth(se = F) + 
  theme_tufte()


with(D, plot(mean, sd))
with(D, plot(mean^2, var))
with(D, plot(sd, intensity))
with(D, plot(var, intensity))

with(D, plot(log(intensity), sd))
with(D, plot(sqrt(intensity), cnt))
with(D, plot(cnt, log(intensity)))


rlm(data=D, log_intensity~cnt) %>% plot

  D %>% 
  filter(cnt >= 5, cnt <=12) %>%
  group_by(cnt) %>%
  summarize(mean_log_intensity = mean(log(intensity))) %>%
  with(points(mean_log_intensity, cnt, col='red'))


WW =D %>% 
    filter(cnt >= 5, cnt <=12) %>%
    group_by(cnt) %>%
    summarize(mean_log_intensity = mean(log(intensity))) 

diff(WW$mean_log_intensity)
  
with(D, plot(log(intensity), log(sd)))
with(D, plot(log(sd), log(intensity), pch='.'))
with(D, plot(log(intensity), sd, pch='.'))
with(D, plot(mean, sd, pch='.'))





with(D %>% filter(log_mean > 6.25,
                  log_sd > -6,
                  cnt > 4, 
                  cnt < 10), 
     plot(log_mean, log_sd, col=cnt))

qplot(x=log_mean, y=log_sd, data=D) + theme_tufte() +
  geom_smooth(method = 'lm', aes(color=ordered(cnt))) 

ggplot(data=D, aes(x=ordered(cnt), y=log_intensity)) + 
  geom_boxplot() + theme_tufte()

ggplot(data=D, aes(x=sd, y=log_intensity)) + 
  geom_point() + theme_tufte()



nice_model = D %>%
  filter(sd > 0) %>%
  rlm(data=., log_sd ~ log_mean + log_intensity + 0) 
plot(nice_model)
summary(nice_model)


D %>% qplot(x = log_intensity, y=log_sd, data=.) + 
  theme_tufte() +
  geom_smooth()



N = lm(sd ~ poly(mean, 2) + poly(intensity, 2), D)

summary(N)
plot(N)
N = rlm(sd ~ poly(mean, 2) + cnt, D)
summary(N)
plot(N)


N = rlm(var ~ intensity + poly(mean, 2) + poly(cnt, 2), D)
summary(N)
plot(N)


with(D, plot(mean, sd))
ablines(N)
