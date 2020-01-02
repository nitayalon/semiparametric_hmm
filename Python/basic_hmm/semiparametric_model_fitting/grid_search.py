import numpy as np
from basic_hmm.semiparametric_model_fitting import baum_welch

class GridSearch:

    def __init__(self, observations):
        self.observations = observations

    def external_grid_search(self):
        p = np.arange(0,1,0.1)
        theta = np.arange(0,1,0.1)
        grid = np.array([(x,y) for x in theta for y in p])
        transition_matrix = np.array([0.95,0.025,0.025,0.1,0.8,0.1,0.025,0.025,0.95]).reshape([3,3])
        llk = []
        for pair in grid:
            bw = baum_welch.RestrictedBaumWelch(pair[0], pair[1], self.observations, transition_matrix)
            gs = bw.restricted_baum_welch()
            llk.append(gs[1])
