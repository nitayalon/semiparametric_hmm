from basic_hmm.mock_hmm import basic_hmm
from basic_hmm.semiparametric_model_fitting import forward
from basic_hmm.semiparametric_model_fitting import backward

hmm = basic_hmm.BasicHMM(3, 0.01, 0.05, 2 / 3, [1 / 3, 1 / 3, 1 / 3], emission_density='laplace', sample=True)

forward_vec = forward.Forward(hmm.observations, hmm.n_states, hmm.transition_matrix, hmm.theta_parameter, hmm.q / 2)
backward_vec = backward.Backward(hmm.observations, hmm.n_states, hmm.transition_matrix, hmm.theta_parameter, hmm.q / 2)
print(forward_vec.compute_normalized_llk())