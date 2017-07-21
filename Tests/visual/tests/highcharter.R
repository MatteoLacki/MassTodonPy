library(highcharter)

mpgman3 <- mpg %>% 
    group_by(manufacturer) %>% 
    summarise(n = n(), unique = length(unique(model))) %>% 
    arrange(-n, -unique) %>% 
    glimpse()
## Observations: 15
## Variables: 3
## $ manufacturer <chr> "dodge", "toyota", "volkswagen", "ford", "chevrol...
## $ n            <int> 37, 34, 27, 25, 19, 18, 14, 14, 13, 9, 8, 5, 4, 4, 3
## $ unique       <int> 4, 6, 4, 4, 4, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1

hchart(mpgman3, "treemap", hcaes(x = manufacturer, value = n, color = unique))


