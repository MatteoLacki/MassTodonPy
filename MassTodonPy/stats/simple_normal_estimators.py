import numpy as np
from math import sqrt

def max_intensity_mean(mz_g, intensity_g):
    """An idiotically simple estimator: the mass of the most intense peak."""
    return mz_g[np.argmax(intensity_g)]


def mean(mz_g, intensity_g):
    """Idiotically simple mean m/z estimator.

    Simply: a weighted mean of m/z."""
    return np.dot(mz_g, intensity_g)/sum(intensity_g)


def sd(mz_g, intensity_g, mean_mz = None):
    """Idiotically simple standard deviation of m/z estimator.

    Simply: a weighted sum of deviations to the mean.."""
    if mean_mz is None:
        mean_mz = mean(mz_g, intensity_g)
    probs   = intensity_g/sum(intensity_g)
    return sqrt(np.dot((mz_g - mean_mz)**2, probs) )
