from ggplot import ggplot, aes, geom_rect, theme_bw, theme_xkcd
from pandas import DataFrame as DF


# X = DF({'mz': range(5), 'intensity':range(5)})
def plot_spectrum(X):
    print ggplot(aes(xmin='mz-.01', xmax='mz+.01', ymin='0', ymax='intensity'), data=X) +\
        geom_rect() +\
        theme_xkcd()
