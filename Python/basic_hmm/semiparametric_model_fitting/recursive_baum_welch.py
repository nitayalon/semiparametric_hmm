from basic_hmm.semiparametric_model_fitting import forward
from basic_hmm.semiparametric_model_fitting import backward
from basic_hmm.semiparametric_model_fitting import compute_central_distribution
import numpy as np


class RecursiveBaumWelch():

    def __init__(self, observations, theta, p, initial_matrix):
        self.observations = observations
        self.n_obs = len(observations)
        self.theta = theta
        self.p = p
        self.initial_matrix = initial_matrix
        self.transition_matrix_dim = initial_matrix.shape[0]
        self.pseudo_count = 0.01

    def compute_polynom_coefficients(self, intermediate_transition_matrix):
        alpha = self.p / (1 - self.p)
        a1 = intermediate_transition_matrix[1, 2] + intermediate_transition_matrix[3, 2] + intermediate_transition_matrix[2, 1] + intermediate_transition_matrix[2, 3]
        a2 = intermediate_transition_matrix[1, 1] + intermediate_transition_matrix[3, 3]
        a3 = intermediate_transition_matrix[2, 2]
        b1 = intermediate_transition_matrix[3, 1] + intermediate_transition_matrix[1, 3]
        b2 = intermediate_transition_matrix[1, 1] + intermediate_transition_matrix[3, 3]
        return [a1, a2, a3, b1, b2, alpha]

    def solve_for_s(self, polynom_coefficients):
        pass

    def baum_welch(self):
        forward_function = forward.Forward(self.observations, self.transition_matrix_dim, self.initial_matrix, self.theta, self.p)
        backward_function = backward.Backward(self.observations, self.transition_matrix_dim, self.initial_matrix, self.theta, self.p)
        central_density = compute_central_distribution.ComputeCentralDistribution(self.observations)
        forward_vector = forward_function.forward_vector()
        backward_vector = backward_function.backward_vector()
        transition_matrix = np.zeros(self.initial_matrix.shape)
        all_probabilities = central_density.compute_all_probabilities(self.theta, self.p)

        prob_observations = forward_vector[1, self.observations.size]
        for i in range(1, self.transition_matrix_dim):
            # last observation, ith row
            j = forward_vector[i, self.observations.size]
            if j > -np.inf:
                prob_observations = j + np.log(1 + np.exp(prob_observations - j))

        for x in range(self.transition_matrix_dim):
            for y in range(self.transition_matrix_dim):
                temp = forward_vector[x, 1] + np.log(self.initial_matrix[x, y]) + \
                   np.log(all_probabilities[y, 1 + 1]) + backward_vector[y, 1 + 1]
                for i in range(1,self.n_obs - 2):
                    j = forward_vector[x, i] + np.log(self.initial_matrix[x, y]) + \
                        np.log(all_probabilities[y, i + 1]) + backward_vector[y, i + 1]
                    if j > -np.inf:
                        temp = j + np.log(1 + np.exp(temp - j))
                temp = np.exp(temp - prob_observations)
                transition_matrix[x, y] = temp

        polynom_coefficients = self.compute_polynom_coefficients(transition_matrix)
        try:
            s = self.solve_for_s(polynom_coefficients)
        except:
            s = 0