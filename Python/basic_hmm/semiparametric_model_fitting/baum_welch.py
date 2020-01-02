import numpy as np
from basic_hmm.semiparametric_model_fitting import recursive_baum_welch
from itertools import count


class RestrictedBaumWelch:

    def __init__(self, theta, p, observations, initial_matrix, number_of_iterations = None):
        self.theta = theta
        self.p = p
        self.observations = observations
        self.transition_matrix = initial_matrix
        self.llk = []
        self.s = []
        self.t = []
        self.delta = 1e-5
        if number_of_iterations is None:
            self.number_of_iterations = number_of_iterations
        else:
            self.number_of_iterations = 5000

    def restricted_baum_welch(self):
        for i in count(0):
            done = False
            rbw = recursive_baum_welch.RecursiveBaumWelch(self.observations, self.theta, self.p, self.transition_matrix)
            bw_results = rbw.baum_welch()
            self.llk.append(-bw_results[1])
            self.s.append(bw_results[2])
            self.t.append(bw_results[3])
            self.transition_matrix = bw_results[0]
            if i > 0:
                done = np.abs(self.llk[i] - self.llk[i-1]) < self.delta
            if done:
                break
            if i > 0 & i % 10 == 0:
                print(f'Iteration number {i}, llk diff = {self.llk[i] - self.llk[i-1]}')
                print(f'Iteration number {i}, s ratio = {self.s[i]/self.s[i-1]}')
        return self.transition_matrix, self.llk[-1]




