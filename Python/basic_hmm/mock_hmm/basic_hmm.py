import numpy as np
from basic_hmm.mock_hmm import generate_samples_from_laplace
from basic_hmm.mock_hmm import generate_samples_from_gaussian


class BasicHMM:

    def __init__(self, n_states, s, t, q, initial_vector,
                 theta_parameter=0.05,
                 emission_density=None):
        if emission_density is None:
            emission_density = 'normal'
        self.n_states = n_states
        self.s = s
        self.t = t
        self.q = q
        self.initial_vector = initial_vector
        self.emission_density = emission_density
        self.theta_parameter = [-theta_parameter, 0.0, theta_parameter]
        self.transition_matrix = self.create_transition_matrix()
        self.hidden_state = []
        self.observations = []

    def create_transition_matrix(self):
        first_row = [np.max(1 - self.s - self.t), self.s, self.t]
        second_row = [self.q / 2, 1 - self.q, self.q / 2]
        third_row = [self.t, self.s, np.max(1 - self.s - self.t)]
        transition_matrix = np.array([first_row, second_row, third_row])
        return transition_matrix

    def generate_hidden_sample(self, sequence_length):
        m = self.transition_matrix.shape[1]
        hidden_states = np.empty(sequence_length)
        hidden_states[0] = np.random.choice(m, 1, p=self.initial_vector)[0]
        for i in range(1, sequence_length):
            hidden_states[i] = np.random.choice(m, 1, p=self.transition_matrix[int(hidden_states[i-1])])
        return hidden_states

    def sample_observations_from_hidden_state(self,hidden_sequence):
        n = len(hidden_sequence)
        observations = np.empty(n)
        for i in range(0, n):
            if self.emission_density == 'laplace':
                obs = generate_samples_from_laplace.sample_from_laplace(1, self.theta_parameter[int(hidden_sequence[i])])
            elif self.emission_density == 'gaussian':
                obs = generate_samples_from_gaussian.sample_from_gaussian(1, self.theta_parameter[int(hidden_sequence[i])])
            else:
                obs = generate_samples_from_gaussian.sample_from_gaussian(1, self.theta_parameter[int(hidden_sequence[i])])
            observations[i] = obs
        return observations

    def generate_observation(self, sequence_length):
        hidden_sequence = self.generate_hidden_sample(sequence_length)
        observation = self.sample_observations_from_hidden_state(hidden_sequence)
        self.hidden_state = hidden_sequence
        self.observations = observation
        return observation

# if __name__ == "__main__":
#     hmm = BasicHMM(3, 0.01, 0.05, 2 / 3, [1 / 3, 1 / 3, 1 / 3], emission_density='laplace')
#     print(hmm.generate_observation(1500))


