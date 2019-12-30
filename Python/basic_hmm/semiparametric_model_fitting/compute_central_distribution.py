import numpy as np


def compute_exponential_kernel(observation, theta, b):
    return np.array([np.exp(-theta * observation - b), 1, np.exp(theta * observation - b)])


class ComputeCentralDistribution:

    def __init__(self, observations):
        self.data = observations
        self.n_obs = observations.size

    def compute_b(self, theta, p):
        b = self.newton_raphson(self.data, theta, p)
        return b

    def newton_raphson(self, data, theta, p, iterations=1000):
        b = np.empty(iterations)
        for i in range(1, iterations):
            b[i] = b[i-1] - self.__optimization_target_function(data, theta, p, b[i-1]) / \
                   self.__compute_newton_raphson_denom(data, theta, p, b[i-1])
        return b[iterations - 1]

    def __optimization_target_function(self, data, theta, p, b):
        n = len(data)
        target_function = self.__compute_central_distribution(data, theta, p, b)
        y = sum(target_function)
        return y-n

    def __compute_central_distribution(self, x, theta, p, b):
        g_x = 1 / (1 - p + p * np.exp(-b) * np.cosh(theta * x))
        return g_x

    def __compute_newton_raphson_denom(self, x, theta, p, b):
        numerator = p * np.exp(-b) * np.cosh(theta * x)
        denom = (1 - p + p * np.exp(-b) * np.cosh(theta * x)) ** 2
        return sum(numerator / denom)

    def exponential_kernel(self, theta, p):
        exponential_kernel = 1 - p + p / 2 * np.exp(-self.compute_b(theta, p)) * np.cosh(theta * self.data)
        return exponential_kernel

    def point_density(self, theta, p):
        point_density = (1 / self.n_obs) * (1 / self.exponential_kernel(theta, p))
        return point_density

    def compute_all_probabilities(self, theta, p):
        central_distribution = np.array(self.point_density(theta, p))
        central_distribution = central_distribution.reshape(10000, 1)
        b = self.compute_b(theta, p)
        exponential_kernel = np.array([compute_exponential_kernel(xi, theta, b) for xi in self.data])
        all_probabilities = np.multiply(central_distribution, exponential_kernel)
        return all_probabilities

