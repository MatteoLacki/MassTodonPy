library(tidyverse)

ubi = read_csv('Ubiquitin_ETD_10 ms_1071.csv') %>% select(mz, intensity)


ubi %>% filter(1420 < mz, mz < 1440) %>% plot(type='h')
ubi %>% filter(1040 < mz, mz < 1080) %>% plot(type='h')

