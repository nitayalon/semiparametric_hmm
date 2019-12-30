import numpy as np
from basic_hmm.semiparametric_model_fitting import compute_central_distribution


class Backward:

    def __init__(self, observations, n_states, transition_matrix, theta, p):
        self.observations = observations
        self.n_states = n_states
        self.transition_matrix = transition_matrix
        self.theta = np.abs(theta[0])
        self.p = p
        self.n_obs = len(observations)
        self.central_distribution = compute_central_distribution.ComputeCentralDistribution(observations)
        self.all_probabilities = self.central_distribution.compute_all_probabilities(self.theta, self.p)

    def backward_vector(self):
        backward_vector = np.zeros(self.all_probabilities.shape)
        for k in range(self.n_obs - 2, -1, -1):
            for state in range(0, self.n_states):
                logsum = -np.inf
                for next_state in range(0, self.n_states):
                    temp = backward_vector[k + 1, next_state] + \
                           np.log(self.transition_matrix[state, next_state] * self.all_probabilities[k + 1, state])
                    if temp > -np.inf:
                        logsum = temp + np.log(1 + np.exp(logsum - temp))
                backward_vector[k, state] = logsum
        return backward_vector

