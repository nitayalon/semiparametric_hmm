from basic_hmm.semiparametric_model_fitting import forward
from basic_hmm.semiparametric_model_fitting import backward
from basic_hmm.semiparametric_model_fitting import compute_central_distribution
import numpy as np


class RecursiveBaumWelch:

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
        a1 = intermediate_transition_matrix[0, 1] + intermediate_transition_matrix[2, 1] + \
             intermediate_transition_matrix[1, 0] + intermediate_transition_matrix[1, 2]
        a2 = intermediate_transition_matrix[0, 0] + intermediate_transition_matrix[2, 2]
        a3 = intermediate_transition_matrix[1, 1]
        a4 = intermediate_transition_matrix[2, 0] + intermediate_transition_matrix[0, 2]
        return [a1, a2, a3, a4,alpha]

    def solve_for_s(self, polynom_coefficients):
        a1, a2, a3, a4, alpha = polynom_coefficients
        upper_limit = np.min([1, 1/alpha])
        ratio = a4 / (a3 + a4)
        A = (alpha - alpha* ratio) * a1 + (ratio - 1)*a2 + alpha * a3
        B = (alpha * (ratio - 1) + (ratio - 1)) * a1 + (1- ratio) * a2 - a3
        C = (1-ratio) * a1
        discriminante = B ** 2 - 4 * A * C
        if discriminante < 0:
            return np.array([0])
        solutions = np.array([[(-B - np.sqrt(discriminante)) / (2 * A)], [(B - np.sqrt(discriminante)) / (2 * A)]])
        if all(solutions < 0):
            return np.array([0])
        if all(solutions > upper_limit):
            return np.array([upper_limit])
        if any(solutions < upper_limit) and any(solutions > 0):
            return self.find_feasible_solution(solutions, upper_limit)
        return np.array([0.5])

    def find_feasible_solution(self, solution_vector, upper_limit):
        if all(solution > 0 for solution in solution_vector) and \
                all(solution < upper_limit for solution in solution_vector):
                return solution_vector
        feasible_value = solution_vector[(solution_vector > 0) & (solution_vector < upper_limit)]
        if len(feasible_value) == 0:
            return upper_limit
        return feasible_value

    def solve_for_t(self, feasible_s_value, polynom_coefficients):
        a1, a2, a3, a4, alpha = polynom_coefficients
        epsilon = 0.01
        a4 = a4 + epsilon
        ratio  = a4 / (a3 + a4)
        current_value = (1 - feasible_s_value) * ratio
        if all(feasible_s_value == 0):
            return ratio
        if all(current_value <= 0):
            return np.array([0])
        if all(current_value >= 1):
            return np.array([1])
        return current_value[(current_value > 0) & (current_value < 1)]

    def find_optimal_transition_matrix(self,s,t):
        s_min, s_max = s
        t_min, t_max = s
        q_min = np.min([self.p / (1 - self.p) * s_min, 1])
        q_max = np.min([self.p / (1 - self.p) * s_max, 1])
        min_list = self.create_transition_matrix(s_min, t_min, q_min)
        max_list = self.create_transition_matrix(s_max, t_max, q_max)
        if min_list[1] > max_list[1]:
            return min_list
        else:
            return max_list

    def create_transition_matrix(self, s, t, q):
        transition_matrix = np.array([1 - s - t, s, t, q / 2, 1 - q, q / 2, t, s, 1 - s - t])
        transition_matrix = transition_matrix.reshape([3, 3])
        transition_matrix = transition_matrix + 1e-4
        transition_matrix = transition_matrix / np.sum(transition_matrix, axis=1)
        forward_vector = forward.Forward(self.observations, self.transition_matrix_dim, transition_matrix, self.theta, self.p)
        llk = forward_vector.compute_normalized_llk()
        return transition_matrix, llk, s, t, forward_vector

    def baum_welch(self):
        forward_function = forward.Forward(self.observations, self.transition_matrix_dim, self.initial_matrix, self.theta, self.p)
        backward_function = backward.Backward(self.observations, self.transition_matrix_dim, self.initial_matrix, self.theta, self.p)
        central_density = compute_central_distribution.ComputeCentralDistribution(self.observations)
        forward_vector = forward_function.forward_vector()
        backward_vector = backward_function.backward_vector()
        transition_matrix = np.zeros(self.initial_matrix.shape)
        all_probabilities = central_density.compute_all_probabilities(self.theta[0], self.p)

        prob_observations = forward_vector[-1, 0]
        for i in range(1, self.transition_matrix_dim):
            # last observation, ith row
            j = forward_vector[-1, i]
            if j > -np.inf:
                prob_observations = j + np.log(1 + np.exp(prob_observations - j))

        for x in range(self.transition_matrix_dim):
            for y in range(self.transition_matrix_dim):
                temp = forward_vector[x, 0] + np.log(self.initial_matrix[x, y]) + \
                   np.log(all_probabilities[y, 1]) + backward_vector[y, 1]
                for i in range(1,self.n_obs - 1):
                    j = forward_vector[i, x] + np.log(self.initial_matrix[x, y]) + \
                        np.log(all_probabilities[i + 1, y]) + backward_vector[i + 1, y]
                    if j > -np.inf:
                        temp = j + np.log(1 + np.exp(temp - j))
                temp = np.exp(temp - prob_observations)
                transition_matrix[x, y] = temp
        polynom_coefficients = self.compute_polynom_coefficients(transition_matrix)
        try:
            s = self.solve_for_s(polynom_coefficients)
        except TypeError:
            s = self.solve_for_s(polynom_coefficients)
            print(f'Something went wrong with s solver')
        t = self.solve_for_t(s, polynom_coefficients)
        if all(s == 0):
            q = self.p / (1 - self.p) * s
            optimal_transition_matrix = self.create_transition_matrix(s,t,q)
        elif len(s) > 1:
            optimal_transition_matrix = self.find_optimal_transition_matrix(s,t)
        else:
            q = self.p / (1 - self.p) * s
            optimal_transition_matrix = self.create_transition_matrix(s,t,q)
        return optimal_transition_matrix[0], optimal_transition_matrix[1], optimal_transition_matrix[2], optimal_transition_matrix[3], optimal_transition_matrix[4]
