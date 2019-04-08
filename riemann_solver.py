import numpy as np

class RiemanSolver:
    def __init__(self, rp):
        
        self.rp = rp
        self.input = rp.input
        self.geo = rp.geo
        self.N = self.geo.N
        self.mat = rp.mat
        self.fields = rp.fields

        self.delta_rho = np.zeros(self.N)
        self.delta_momentum = np.zeros(self.N)
        self.delta_mat_energy = np.zeros(self.N)

        self.rho_L = np.zeros(self.N)
        self.rho_L_new = np.zeros(self.N)

        self.rho_R = np.zeros(self.N)
        self.rho_R_new = np.zeros(self.N)

        self.momentum_L = np.zeros(self.N)
        self.momentum_L_new = np.zeros(self.N)

        self.momentum_R = np.zeros(self.N)
        self.momentum_R_new = np.zeros(self.N)

        self.mat_energy_L = np.zeros(self.N)
        self.mat_energy_L_new = np.zeros(self.N)

        self.mat_energy_R = np.zeros(self.N)
        self.mat_energy_R_new = np.zeros(self.N)

    def computeSlopes(self):

        w = self.input.omega

        rho_old = self.fields.rho_old
        rho_bL = self.fields.rho_bL
        rho_bR = self.fields.rho_bR

        momentum_old = self.fields.momentum_old
        momentum_bL = self.fields.momentum_bL
        momentum_bR = self.fields.momentum_bR

        mat_energy_old = self.fields.mat_energy_old
        mat_energy_bL = self.fields.mat_energy_bL
        mat_energy_bR = self.fields.mat_energy_bR

        self.delta_rho[0] = 0.5 * (1 - w) * (rho_old[0] - rho_bL) + \
                            0.5 * (1 - w) * (rho_old[1] - rho_old[0])

        self.delta_momentum[0] = 0.5 * (1 - w) * (momentum_old[0] - momentum_bL) + \
                                 0.5 * (1 - w) * (momentum_old[1] - momentum_old[0])

        self.delta_mat_energy[0] = 0.5 * (1 - w) * (mat_energy_old[0] - mat_energy_bL) + \
                                   0.5 * (1 - w) * (mat_energy_old[1] - mat_energy_old[0])

        for i in range(1, self.N - 1):

            self.delta_rho[i] = 0.5 * (1 - w) * (rho_old[i] - rho_old[i-1]) + \
                                0.5 * (1 - w) * (rho_old[i+1] - rho_old[i])

            self.delta_momentum[i] = 0.5 * (1 - w) * (momentum_old[i] - momentum_old[i-1]) + \
                                     0.5 * (1 - w) * (momentum_old[i+1] - momentum_old[i])

            self.delta_mat_energy[i] = 0.5 * (1 - w) * (mat_energy_old[i] - mat_energy_old[i-1]) + \
                                       0.5 * (1 - w) * (mat_energy_old[i+1] - mat_energy_old[i])

        self.delta_rho[-1] = 0.5 * (1 - w) * (rho_old[len(rho_old) - 1] - rho_old[len(rho_old) - 2]) + \
                             0.5 * (1 - w) * (rho_bR - rho_old[len(rho_old) - 1])

        self.delta_momentum[-1] = 0.5 * (1 - w) * (momentum_old[len(momentum_old) - 1] - momentum_old[len(momentum_old) - 2]) + \
                             0.5 * (1 - w) * (rho_bR - momentum_old[len(momentum_old) - 1])

        self.delta_mat_energy[-1] = 0.5 * (1 - w) * (mat_energy_old[len(mat_energy_old) - 1] - mat_energy_old[len(mat_energy_old) - 2]) + \
                             0.5 * (1 - w) * (rho_bR - mat_energy_old[len(mat_energy_old) - 1])

    def computeEdgeValues(self):

        rho_old = self.fields.rho_old
        momentum_old = self.fields.momentum_old
        mat_energy_old = self.fields.mat_energy_old

        for i in range(self.N):

            self.rho_L[i] = rho_old[i] - 0.5 * self.delta_rho[i]
            self.rho_R[i] = rho_old[i] + 0.5 * self.delta_rho[i]

            self.momentum_L[i] = momentum_old[i] - 0.5 * self.delta_momentum[i]
            self.momentum_R[i] = momentum_old[i] + 0.5 * self.delta_momentum[i]

            self.mat_energy_L[i] = mat_energy_old[i] - 0.5 * self.delta_mat_energy[i]
            self.mat_energy_R[i] = mat_energy_old[i] + 0.5 * self.delta_mat_energy[i]

    def computeNewEdgeValues(self, dt):

        dr = self.geo.dr

        rho_L = self.rho_L
        rho_R = self.rho_R

        momentum_L = self.momentum_left
        momentum_R = self. momentum_right

        mat_energy_L = self.mat_energy_L
        mat_energy_R = self.mat_energy_R

        for i in range(self.N):

            F_rho_L = momentum_L[i]
            F_rho_R = momentum_R[i]

            P_L = (self.mat.gamma - 1) * rho_L[i] * (mat_energy_L[i] / rho_L[i] - 0.5 * (momentum_L[i] / rho_L[i])**2)
            P_R = (self.mat.gamma - 1) * rho_R[i] * (mat_energy_R[i] / rho_R[i] - 0.5 * (momentum_R[i] / rho_R[i])**2)

            F_momentum_L = momentum_L[i]**2 / rho_L[i] + P_L
            F_momentum_R = momentum_R[i]**2 / rho_R[i] + P_R

            F_mat_energy_L = (momentum_L[i] / rho_L[i]) * (mat_energy_L[i] + P_L)
            F_mat_energy_R = (momentum_R[i] / rho_R[i]) * (mat_energy_R[i] + P_R)

            rho_L_new[i] = rho_L + 0.5 * dt / dr[i] * (F_rho_L - F_rho_R)
            rho_R_new[i] = rho_R + 0.5 * dt / dr[i] * (F_rho_L - F_rho_R)

            momentum_L_new[i] = momentum_L + 0.5 * dt / dr[i] * (F_momentum_L - F_momentum_R)
            momentum_R_new[i] = momentum_R + 0.5 * dt / dr[i] * (F_momentum_L - F_momentum_R)

            mat_energy_L_new[i] = mat_energy_L + 0.5 * dt / dr[i] * (F_mat_energy_L - F_mat_energy_R)
            mat_energy_R_new[i] = mat_energy_L + 0.5 * dt / dr[i] * (F_mat_energy_L - F_mat_energy_R)



        

