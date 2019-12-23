import numpy as np


def sample_from_laplace(sample_size, location = 0):
    if np.abs(location) > 1:
        location = 0.99
    ind = np.random.uniform(0.0, 1.0, sample_size)
    if location < 0:
        direction = -1
    else:
        direction = 1
    theta = np.abs(location)
    y = np.random.exponential(1,sample_size)
    h_x = 1 / (1 + theta) * y
    g_x = 1 / (1 - theta) * y
    p_positive = (1 - theta) / 2
    positive_ind = ind > p_positive
    if direction == -1:
        negative_sign = [1, -1]
    else:
        negative_sign = [-1, 1]
    x = negative_sign[0] * g_x * positive_ind + negative_sign[1] * h_x * (1 - positive_ind)
    return x

