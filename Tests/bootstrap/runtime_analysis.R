library(tidyverse)

load(file = 'runtime_data_ubi.rda')
load(file = 'runtime_data_subP.rda')


runtime_data = bind_rows(
    runtime_data_subP %>% mutate(Q = 3, molecule = 'substance P'),
    runtime_data_ubi  %>% mutate(Q = round(8500/precursorMZ), molecule = 'ubiquitin') 
)

base_breaks <- function(n = 10){
    function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

runtime_plot = 
    runtime_data %>%
    ggplot() +
    geom_boxplot(
        aes( x = factor(Q), y = total_time, color = real_or_bootstrap)
    ) +
    theme_minimal() +
    # theme_solarized_2(light = FALSE) +
    xlab("Initial Charge") +
    ylab('Runtime [sec]') +
    scale_y_continuous(
        trans  = scales::log10_trans(),
        breaks = base_breaks(),
        labels = prettyNum
    ) +
    theme(legend.position="top") + 
    coord_flip() +
    scale_color_ptol(
        labels = c('Bootstrap', 'Real Data') 
    ) +
    guides(color=guide_legend(title=NULL))


cowplot::ggsave(
    '/Users/matteo/Documents/MassTodon/Paper/images/runtime_plot.pdf',
    runtime_plot,
    limitsize = F,
    dpi=200,
    width = 100, height = 80, units = 'mm', scale=1)



