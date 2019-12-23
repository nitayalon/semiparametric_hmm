import numpy as np


def sample_from_gaussian(sample_size, location=0):
    x = np.random.normal(location, 1.0, sample_size)
    return x

