# install.packages("rbokeh")
library(rbokeh)

 %>% ly_points(1:10)

rectangles = tibble( xleft = 1:10, ybottom = 1:10, xright = 2:11, ytop = (2:11)/2)

figure( data = rectangles ) %>%
ly_rect(xleft = xleft*2, ybottom = ybottom, xright = xright, ytop = ytop, color='red') 



figure() %>%
ly_rect(
    xleft, ybottom, xright, ytop, data = figure_data(fig), color = NULL, alpha = 1, hover = NULL, url = NULL, legend = NULL, lname = NULL, lgroup = NULL, ...)