library('tidyverse'); library('stringr'); library('car')


D = list.files('boot_subP_24_07_2017_csv', full.names = T) %>%
    lapply(read.csv) %>%
    lapply(as_data_frame) %>%
    bind_rows() %>%
    mutate(WH = factor(WH), WV = factor(WV) )


D = D %>% filter(ID != 13)
D150 = D %>% filter(WH == '150') %>%
    select(
        ID, WH, WV, 
        basic.prob.ETnoD, real_or_bootstrap,
        basic.count.ETnoD, basic.count.PTR,
        basic.prob.anion_approached_cation,
        basic.count.unreacted_precursors, 
        basic.count.total_reactions
    ) %>%
    split(.$real_or_bootstrap)


capFirst = function(s) 
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")

# Testing if missing values can be attributed to lack of estimates of ETnoD and PTR.
test = D150$boot %>% filter( is.na(basic.prob.ETnoD) ) 
assertthat::are_equal(all(test$basic.count.ETnoD == 0) & all(test$basic.count.PTR == 0), T)


estimates_probs_plot150 =
    ggplot() +
    theme_minimal() +
    geom_boxplot(
        aes(x=WV, y=basic.prob.ETnoD),
        data = D150$boot
    ) +
    geom_point(
        aes(x=WV, y=basic.prob.ETnoD),
        data  = D150$real, 
        color = 'red', 
        shape = 19 
    ) +
    theme( axis.text.x = element_text(angle = 90, hjust = 1) ) +
    scale_y_continuous(
        labels   = scales::percent,
        sec.axis = sec_axis(
            trans = ~1-.,
            labels=scales::percent,
            name  = 'PTR'
        )
    )+
    xlab('Wave Velocity (Wave Height set to 1,5)') +
    ylab('ETnoD')


D150 = D150 %>% lapply(
    function(x)
        x %>% select(
            ID, WH, WV,  
            basic.count.ETnoD, basic.count.PTR
        ) %>%
        gather('reaction', 'count', 4:5)
)


MAX_COUNT = D150 %>% sapply(function(x) max(x$count)) %>% max
BREAKS = seq(0, 1400, by=200)^2
width = .75

estimates_counts_plot150 =
    ggplot() +
    theme_minimal() +
    theme( axis.text.x = element_text(angle = 90, hjust = 1) ) +
    geom_boxplot( 
        aes(x=WV, y=count, fill=reaction, color=reaction),
        data = D150$boot,
        position = position_dodge(width = width)
    ) +
    geom_line( 
        aes(x=WV, y=count, group=reaction, color=reaction),
        data = D150$real,
        alpha= .5,
        position = position_dodge(width = width)
    ) +
    scale_y_continuous(
        trans  = scales::sqrt_trans(),
        breaks = BREAKS,
        limits = c(0, MAX_COUNT)
    )+
    scale_fill_discrete(
        labels = c('ETnoD','PTR')
    ) +
    scale_color_discrete(
        labels = c('ETnoD','PTR')
    ) +
    xlab('Wave Velocity (Wave Height set to 1,5)') +
    ylab('Intensity [square root scale]') +
    theme(legend.position="bottom") 
    
estimates_counts_plot150_no_guides = 
    estimates_counts_plot150 + guides(color=F, fill=F)

################################################################################################


D300 = D %>% 
    filter(WV == '300') %>%
    select( 
        ID, WH, WV, 
        basic.prob.ETnoD, real_or_bootstrap,
        basic.count.ETnoD, basic.count.PTR,
        basic.prob.anion_approached_cation,
        basic.count.unreacted_precursors, 
        basic.count.total_reactions
    ) %>%    
    split(.$real_or_bootstrap)


estimates_probs_plot300 =
    ggplot() +
    theme_minimal() +
    geom_boxplot(
        aes(x=WH, y=basic.prob.ETnoD),
        data = D300$boot
    ) +
    geom_point(
        aes(x=WH, y=basic.prob.ETnoD),
        data  = D300$real, 
        color = 'red', 
        shape = 19 
    ) +
    theme( axis.text.x = element_text(angle = 90, hjust = 1) ) +
    scale_y_continuous(
        labels   = scales::percent,
        sec.axis = sec_axis(
            trans = ~1-.,
            labels=scales::percent,
            name  = 'PTR'
        )
    )+
    xlab('Wave Height (Wave Velocity set to 300)') +
    ylab('ETnoD')

D300 = D300 %>% lapply(
    function(x)
        x %>% select(
            ID, WH, WV,  
            basic.count.ETnoD, basic.count.PTR
        ) %>%
        gather('reaction', 'count', 4:5)
)

width = .75
estimates_counts_plot300 =
    ggplot() +
    theme_minimal() +
    theme( 
        axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    geom_boxplot( 
        aes(x=WH, y=count, fill=reaction, color=reaction),
        data = D300$boot,
        position = position_dodge(width = width)
    ) +
    geom_line( 
        aes(x=WH, y=count, group=reaction, color=reaction),
        data = D300$real,
        alpha= .5,
        position = position_dodge(width = width)
    ) +
    scale_fill_discrete(
        labels = c('ETnoD','PTR')
    ) +
    scale_color_discrete(
        labels = c('ETnoD','PTR')
    ) +
    scale_y_continuous(
        trans  = scales::sqrt_trans(),
        breaks = BREAKS,
        limits = c(0, MAX_COUNT),
        position = "right"
    )+
    xlab('Wave Height (Wave Velocity set to 300)') +
    ylab('Intensity [square root scale]') 


estimates_counts_plot300_no_guides = 
    estimates_counts_plot300 +
    guides(fill=F, color=F)

# WH_WV = cowplot::plot_grid(
#     estimates_counts_plot150, 
#     estimates_counts_plot300, 
#     ncol=2, rel_widths = c(1.5,.6))


legend <- cowplot::get_legend(estimates_counts_plot150)

WH_WV = cowplot::plot_grid(
    estimates_counts_plot150_no_guides, 
    estimates_counts_plot300_no_guides, 
    ncol=2, align='h', rel_widths = c(1.5,.6))

final_plot = cowplot::plot_grid(WH_WV, legend, nrow=2, rel_heights = c(1,.1))


cowplot::ggsave(
    'substance_P_intensities_of_etnod_ptr.pdf',
    final_plot, limitsize = F, dpi=100, width = 400, height = 100, units = 'mm', scale=1)


cowplot::ggsave(
    'substance_P_intensities_of_etnod_ptr.png',
    final_plot, limitsize = F, dpi=100, width = 400, height = 100, units = 'mm', scale=1)

cowplot::ggsave(
    'substance_P_intensities_of_etnod_ptr_300.png',
    final_plot, limitsize = F, dpi=300, width = 400, height = 100, units = 'mm', scale=1)

#### runtime analysis

D %>% colnames

runtime_data_subP =
    D %>% 
    select( ID, WH, WV, real_or_bootstrap,
            total_time )

save(runtime_data_subP, file = 'runtime_data_subP.rda')
