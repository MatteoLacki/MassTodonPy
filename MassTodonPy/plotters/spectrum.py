import matplotlib.pyplot as plt
import numpy as np


def plot_spectrum(mz, intensity,
                  clusters  = None,
                  plt_style = 'dark_background',
                  show      = True):
    """Make a simple visualization of the data.

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
    show : logical
        Immediately show the plot? Alternatively, just add it to the present canvas.

    """
    plt.vlines(x=mz, ymin=[0], ymax=intensity)
    if clusters:
        plt.scatter(x=mz, y=np.zeros(len(mz)), c=clusters)
    if show:
        plt.show()