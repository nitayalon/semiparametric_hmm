import numpy as np
from basic_hmm.mock_hmm import basic_hmm

class SemiParametricModel:

    def __init__(self, observations):
        self.data = observations

    def compute_b(self, theta, p):
        b = self.newton_raphson(self.data, theta, p)
        return b

    def newton_raphson(self, data, theta, p, iterations=1000):
        b = np.empty(iterations)
        for i in range(1, iterations):
            b[i] = b[i-1] - self._optimization_target_function(data, theta, p, b[i-1]) / \
                   self._compute_newton_raphson_denom(data, theta, p, b[i-1])
        return b[iterations]

    def _optimization_target_function(self, data, theta, p, b):
        n = len(data)
        target_function = self._compute_central_distribution(data, theta, p, b)
        y = sum(target_function)
        return y-n

    def _compute_central_distribution(self, x, theta, p, b):
        g_x = 1 / (1 - p + p * np.exp(-b) * np.cosh(theta * x))
        return g_x

    def _compute_newton_raphson_denom(self, x, theta, p, b):
        numerator = p * np.exp(-b) * np.cosh(theta * x)
        denom = (1 - p + p * np.exp(-b) * np.cosh(theta * x)) ** 2
        return sum(numerator / denom)

if __name__ == "__main__":
    hmm = basic_hmm.BasicHMM(3, 0.01, 0.05, 2 / 3, [1 / 3, 1 / 3, 1 / 3], emission_density='laplace')
    obs = hmm.generate_observation(1500)
    spm = SemiParametricModel(obs)
    spm.compute_b(0.05, 1/3)
