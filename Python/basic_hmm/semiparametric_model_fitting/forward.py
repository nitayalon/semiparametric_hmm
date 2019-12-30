import numpy as np
from basic_hmm.semiparametric_model_fitting import compute_central_distribution


class Forward:

    def __init__(self, observations, n_states, transition_matrix, theta, p):
        self.observations = observations
        self.n_states = n_states
        self.transition_matrix = transition_matrix
        self.theta = np.abs(theta[0])
        self.p = p
        self.n_obs = len(observations)
        self.central_distribution = compute_central_distribution.ComputeCentralDistribution(observations)
        self.all_probabilities = self.central_distribution.compute_all_probabilities(self.theta, self.p)

    def forward_vector(self):
        forward_vector = np.zeros(self.all_probabilities.shape)
        forward_vector[0, :] = np.log(np.array([1/3, 1/3, 1/3]) * self.all_probabilities[0, :])
        for k in range(1, self.n_obs):
            for state in range(0, self.n_states):
                logsum = -np.inf
                for previous_state in range(0, self.n_states):
                    temp = forward_vector[(k - 1), previous_state] + \
                           np.log(self.transition_matrix[previous_state, state])
                    if temp > -np.inf:
                        logsum = temp + np.log(1 + np.exp(logsum - temp))
                forward_vector[k, state] = np.log(self.all_probabilities[k, state]) + logsum
        return forward_vector

