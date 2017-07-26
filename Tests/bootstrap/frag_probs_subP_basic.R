library('tidyverse'); library('stringr'); library('car'); library(ggthemes)

D = list.files('boot_subP_24_07_2017_csv', full.names = T) %>%
    lapply(read.csv) %>%
    lapply(as_data_frame) %>%
    bind_rows() %>%
    mutate(WH = factor(WH), WV = factor(WV) )

D = D %>% filter(ID != 13)

D150 = D %>% filter(WH == '150')
D300 = D %>% filter(WV == '300')

fasta = strsplit('RPKPQQFFGLM', '')[[1]]
get_AA = function(x)
    fasta[x]

uniform_prob = 1/sum(fasta != "P")

# D_no = D150
prep_data = function(D_no) D_no %>%
    select( WH, WV, ID, real_or_bootstrap, 
            matches('basic.prob.'), 
    ) %>%
    select(-contains('ETnoD')) %>%
    select(-contains('PTR')) %>%
    select(-contains('fragmentation')) %>%
    select(-contains('anion')) %>%
    gather( 'AA', 'prob', -c(1:4), na.rm = TRUE ) %>%
    mutate(
        AA_no= as.integer(str_match(AA, "[0-9]+"))+1,
        AA   = get_AA(AA_no)
    ) %>% split(.$real_or_bootstrap)

D4P150 = prep_data(D150)

names(fasta) = 1:length(fasta)
to_string    = as_labeller(fasta)

plot150 =
    ggplot( 
        data = D4P150$boot,
        aes(x=WV, y=prob)
    ) +
    geom_hline(
        aes( yintercept = uniform_prob), 
        linetype = "dashed",
        color    = "orange"
    )+
    geom_boxplot() +
    theme_minimal()+
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(
        labels   = scales::percent
    ) +
    theme(
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    facet_grid(
        .~factor(AA_no, levels=1:length(fasta)),
        drop     = FALSE,
        labeller = to_string
    ) +
    xlab('Wave Velocity (WH = 1,5)') +
    ylab('Fragmentation Probability') +
    coord_flip()



###########################

D4P300 = prep_data(D300)
plot300=
    ggplot( 
        data = D4P300$boot,
        aes(x=WH, y=prob)
    ) +
    geom_hline(
        aes(yintercept=uniform_prob), 
        linetype = "dashed",
        color    = "orange"
    )+
    geom_boxplot() +
    theme_minimal() +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(
        labels = scales::percent
    ) +
    theme(
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    facet_grid(
        .~factor(AA_no, levels=1:length(fasta)),
        drop     = FALSE,
        labeller = to_string
    ) +
    xlab('Wave Height (WV = 300)') +
    ylab('Fragmentation Probability') +
    coord_flip()+
    guides(fill=FALSE, color=FALSE)


############################ combine plots

big_plot = cowplot::plot_grid(plot150, plot300, nrow=2, align='v', rel_heights = c(1.5,.65) )

date = format(Sys.time(), "%a_%b_%d_%X_%Y")
cowplot::ggsave(
    paste0('/Users/matteo/Documents/MassTodon/Paper/images/Fragmentation_plot-', date, '.pdf', sep='', collapse=''),
    big_plot,
    limitsize = F,
    dpi=600,
    width = 300, height = 300, units = 'mm', scale=1)

date = format(Sys.time(), "%a_%b_%d_%X_%Y")
cowplot::ggsave(
    paste0('/Users/matteo/Documents/MassTodon/Paper/images/Fragmentation_plot-', date, '.png', sep='', collapse=''),
    big_plot,
    limitsize = F,
    dpi=200,
    width = 300, height = 300, units = 'mm', scale=1)

# ggplot( data = D4P$boot %>% rename(Algorithm=algo),
#         aes(x=WV, y=prob, fill=Algorithm, color=Algorithm)
# ) +
#   geom_hline(aes(yintercept=uniform_prob), linetype="dashed")+
#   geom_boxplot() +
#   theme_minimal() +
#   scale_x_discrete(drop=FALSE) +
#   scale_y_continuous(
#     labels   = scales::percent
#   ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(.~factor(AA_no, levels=1:length(fasta)), drop=FALSE)+
#   xlab('Wave Velocity (Wave Height set to 1,5)') +
#   ylab('Fragmentation Probability') +
#   coord_flip()

# ggplot( data = D4P$boot,
#         aes(x=WV, y=prob, fill=algo, color=algo)
# ) +
# geom_boxplot() +
# theme_minimal() +
# scale_x_discrete(drop=FALSE) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# facet_grid(AA_no~.)
#
#
#
# ggplot( data = D4P$boot,
#         aes(x=WV, y=prob, color=algo, fill=algo)
# ) +
#   geom_boxplot() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(factor(AA_no, levels=1:length(fasta))~.)
#
#
#
# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob, color=algo, fill=algo)
# ) +
#   geom_boxplot() +
#   scale_x_discrete(drop=FALSE) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   coord_flip()+
#   facet_grid(WV~.)
#
#
# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob)
# ) +
#   geom_boxplot() +
#   scale_x_discrete(drop=FALSE) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(.~WV)
#


# ggplot( data = D4P$boot,
#         aes(x=WV, y=prob, fill=algo, color=algo)
# ) +
#   geom_boxplot() +
#   theme_minimal() +
#   scale_x_discrete(drop=FALSE) +
#   scale_y_continuous(
#     labels   = scales::percent
#   ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(factor(AA_no, levels=1:length(fasta))~., drop=FALSE)+
#   xlab('Wave Velocity (Wave Height set to 1,5)') +
#   ylab('Fragmentation Probability')



# binary search in the space of all fucking plots...

# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob, fill=factor(WV), color=factor(WV))
# ) +
# geom_boxplot() +
# theme_minimal() +
# scale_x_discrete(drop=FALSE) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# facet_grid(algo~.)
#
#
# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob)
# ) +
# geom_boxplot() +
# theme_minimal() +
# scale_x_discrete(drop=FALSE) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# facet_grid(algo~WV)

# Looks like a dirty toilet
# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob)
# ) +
#   geom_boxplot() +
#   theme_minimal() +
#   scale_x_discrete(drop=FALSE) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(WV~algo, scales="free_y")





# ggplot( data = D4P$boot,
#         aes(x=algo, y=prob)
# ) +
# geom_boxplot() +
# theme_minimal() +
# scale_x_discrete(drop=FALSE) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# facet_grid(WV~AA_no)


#
# ggplot( data = D4P$boot,
#         aes(x=factor(AA_no, levels=1:length(fasta)), y=prob), fill=algo, color=algo
# ) +
# geom_boxplot() +
# theme_minimal() +
# scale_x_discrete(drop=FALSE) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# facet_grid(WV~., scales="free_y")
#
