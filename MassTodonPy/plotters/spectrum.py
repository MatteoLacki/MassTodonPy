import matplotlib.pyplot as plt
import numpy as np


def plot_spectrum(mz, intensity,
                  clusters  = None,
                  plt_style = 'dark_background',
                  colors_no = 10,
                  show      = True):
    """Make a simple visualization of a mass spectrum.

    If provided, the clusters will be represented as 
    colorful dots on the bottom of the peaks.

    Parameters
    ----------
    mz : np.array
        The recorded m/z values of the spectrum.
    intensity : np.array
        The recorded intensities of the spectrum.
    clusters : np.array
        The assignments into different peak clusters.
    plt_style : str
        The type of the matplotlib style used in the plot.
    colors_no : int
        The number of colors to cycle through.
    show : logical
        Immediately show the plot? Alternatively, just add it to the present canvas.

    """
    plt.style.use(plt_style)
    plt.vlines(x=mz, ymin=[0], ymax=intensity)
    if clusters is not None:
        cmap = plt.get_cmap('tab10', colors_no)
        colors = [cmap(c % colors_no) for c in clusters]
        plt.scatter(x=mz, y=np.zeros(len(mz)), c=colors)
    if show:
        plt.show()