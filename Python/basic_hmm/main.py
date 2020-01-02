from basic_hmm.mock_hmm import basic_hmm
from basic_hmm.semiparametric_model_fitting import baum_welch
from basic_hmm.semiparametric_model_fitting import forward


if __name__ == '__main__':
    hmm = basic_hmm.BasicHMM(3, 0.01, 0.05, 2 / 3, [1 / 3, 1 / 3, 1 / 3], emission_density='normal', sample=True)
    forward_vec = forward.Forward(hmm.observations, 3, hmm.transition_matrix, hmm.theta_parameter, hmm.q / 2)
    original_llk = forward_vec.compute_normalized_llk()
    bw = baum_welch.RestrictedBaumWelch(hmm.theta_parameter, hmm.q / 2, hmm.observations,  hmm.transition_matrix)
    gs = bw.restricted_baum_welch()
    print(gs[1], -original_llk)
    print(gs[0])
