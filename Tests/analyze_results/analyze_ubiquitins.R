library(tidyverse)

initial_run_files = list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun') %>% tools::file_path_sans_ext()
D = read.csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_only_real/raw_estimates.csv') %>% tbl_df


ubi_overview =
    D %>% 
    filter(molType == 'precursor') %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point() +
    facet_grid( retentionTime ~ interaction(supActive, preActive) )


ubi_overview_2 =
    D %>% 
    filter(molType == 'precursor') %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point() +
    facet_wrap( ~files )



D %>% 
    filter(molType == 'precursor') %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) %>%
    filter( files == 'Ubiquitin_ETD_952_supp_act_70 ms' ) %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point() 



cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/ubi_overview.pdf', ubi_overview, limitsize = F, dpi=100, width = 250, height = 400, units = 'mm', scale=1)
cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/ubi_overview_2.pdf', ubi_overview_2, limitsize = F, dpi=100, width = 250, height = 400, units = 'mm', scale=1)




ubi_6_data = 
    D %>% 
    filter(molType == 'precursor') %>%
    filter(precursor_charge == 6) %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) 

ubi_overview_only_Q_6 = ubi_6_data %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point() +
    facet_wrap( ~files )


W = 
ubi_6_data %>% group_by(files) %>% 
summarize(
    ETnoD_cnt = sum(estimate*ETnoD),
    PTR_cnt   = sum(estimate*PTR)
) %>%
mutate(
    ETnoD_prob = ETnoD_cnt/(ETnoD_cnt + PTR_cnt),
    ETnoD_prob = PTR_cnt/(ETnoD_cnt + PTR_cnt)
) 

W$ETnoD_prob %>% plot







ubi_985_data = 
    D %>% 
    filter(molType == 'precursor') %>%
    filter(precursor_charge == 9) %>%
    filter( files == 'Ubiquitin_ETD_952_supp_act_70 ms.mzXML' ) %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) 

ubi_overview_only_Q_6 = ubi_985_data %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point() +
    facet_wrap( ~files )


W = ubi_985_data %>%
    summarize(
        ETnoD_cnt = sum(estimate*ETnoD),
        PTR_cnt   = sum(estimate*PTR)
    ) %>%
    mutate(
        ETnoD_prob = ETnoD_cnt/(ETnoD_cnt + PTR_cnt),
        PTR_prob = PTR_cnt/(ETnoD_cnt + PTR_cnt)
    ) 





ubi_1420_data = 
    D %>% 
    filter(molType == 'precursor') %>%
    filter(precursor_charge == 6) %>%
    filter( files == 'Ubiquitin_ETD_1428_2 ms.mzXML' ) %>%
    mutate( ETnoD = g, PTR = precursor_charge - g - q, files = tools::file_path_sans_ext(files) ) %>%
    filter( !(files %in% initial_run_files) ) 


ubi_overview_only_Q_6 = 
    ubi_1420_data %>%
    ggplot(aes(x = PTR, y = ETnoD, size = estimate)) + 
    theme_minimal()+
    geom_point()


W = ubi_1420_data %>%
    summarize(
        ETnoD_cnt = sum(estimate*ETnoD),
        PTR_cnt   = sum(estimate*PTR)
    ) %>%
    mutate(
        ETnoD_prob = ETnoD_cnt/(ETnoD_cnt + PTR_cnt),
        PTR_prob = PTR_cnt/(ETnoD_cnt + PTR_cnt)
    ) 

ubispec = read.csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/dupa.csv') %>% tbl_df() %>%
    rename(mz = X0, intensity = X )
ubispec = ubispec[,2:1]
ubispec %>% plot(type='h')
ubispec %>% filter((8563.63 + 6 )/ 6 - 4 < mz,  (8563.63 + 6 )/ 6 + 4 > mz) %>% plot(type='h')

real = ubispec %>% filter(mz < 1431, mz > 1427.5) 
real = ubispec %>% filter(mz < 1716, mz > 1713)
real = ubispec %>% filter((8563.63 + 6 )/ 4 - 4 < mz,  (8563.63 + 6 )/ 4 + 4 > mz) %>% plot(type='h')
real = ubispec %>% filter(mz < 2144, mz > 2141)
real$tag = 'real'

library(IsoSpecR)
ptr = IsoSpecify(molecule = c(C=378, H=629, N=105,O=118,S=1), stopCondition = .99 ) %>% data.frame %>% mutate( intensity = exp(logProb) ) %>% 
    mutate(mz = (mass + 5)/4 ) %>% 
    mutate( mz = round(mz, 1) )%>% 
    group_by( mz ) %>%
    summarize( intensity = sum(intensity) ) %>%
    select( mz, intensity )
ptr$tag = 'ptr'

etnod = IsoSpecify(molecule = c(C=378, H=629, N=105,O=118,S=1), stopCondition = .99 ) %>% data.frame %>% mutate( intensity = exp(logProb) ) %>% 
    mutate(mz = (mass + 6)/4 ) %>% 
    mutate( mz = round(mz, 1) )%>% 
    group_by( mz ) %>%
    summarize( intensity = sum(intensity) ) %>%
    select( mz, intensity )
etnod$tag = 'etnod'

spectra = bind_rows(real, ptr, etnod)
# spectra = bind_rows(real, etnod)
spectra %>% group_by(tag) %>% mutate(intensity = intensity/sum(intensity)) %>% ungroup %>%
    ggplot(aes(x = mz, y = intensity, color = tag )) +
    geom_point()

spectra %>% filter( mz > 2142.8, mz < 2143) 

