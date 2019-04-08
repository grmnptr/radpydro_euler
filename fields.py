import numpy as np

class Fields:
    def __init__(self, rp):
        self.rp = rp
        self.input = rp.input
        self.geo = rp.geo
        self.mat = rp.mat
        self.N = rp.geo.N

        # Densities (piece-wise constant)
        self.rho = self.initializeAtCenters(self.input.rho)
        self.rho_old = np.copy(self.rho)

        # Velocities (piece-wise constant)
        self.u = self.initializeAtCenters(self.input.u)
        self.u_old = np.copy(self.u)

        # Set velocity values in the ghost cells
        self.u_bL = self.input.u_bL
        self.u_bR = self.input.u_bR

        # Internal energy (piece-wise constant)
        self.e = self.initializeAtCenters(self.input.T)
        self.e_old = np.copy(self.e)

        # Set internal energy values in the ghost cells
        self.e_bL = self.input.e_bL
        self.e_bR = self.input.e_bR

        # Pressures (piece-wise constant)
        self.P = (self.mat.gamma - 1) * self.e * self.rho
        self.P_old = np.copy(self.P)

        # Set pressure values in the ghost cells
        self.P_bL = (self.mat.gamma - 1) * self.e_bL * self.rho_bL
        self.P_bR = (self.mat.gamma - 1) * self.e_bR * self.rho_bR

        # Temperatures (piece-wise constant)
        self.T = self.e / self.mat.C_v
        self.T_old = np.copy(self.T)

        # Set pressure values in the ghost cells
        self.T_bL = self.e_bL / self.mat.C_v
        self.T_bR = self.e_bR / self.mat.C_v

        # Compute initial momentum values
        self.momentum = self.rho * self.
        self.momentum_old = np.copy(self.momentum)

        # Set the valuesin ghost cells
        self.momentum_bL = self.rho_bL * self.u_bL
        self.momentum_bR = self.rho_bR * self.u_bR

        # Compute the intiial value of the momentum
        self.mat_energy = self.rho * (0.5 * self.u**2 + self.e)
        self.mat_energy_old = np.copy(self.mat_energy)

        # Setting the values in ghost cells
        self.mat_energy_bR = self.rho_bR * (0.5 * self.u_bR**2 + self.e_bR)
        self.mat_energy_bL = self.rho_bL * (0.5 * self.u_bL**2 + self.e_bL)

        # Fluxes
        self.F_rho = np.zeros(self.N + 1)
        self.F_momentum = np.zeros(self.N + 1)
        self.F_mat_energy = np.zeros(self.N + 1)

    # Copy over all new fields to old positions
    def stepFields(self):
        np.copyto(self.u_old, self.u)
        np.copyto(self.T_old, self.T)
        np.copyto(self.rho_old, self.rho)
        np.copyto(self.P_old, self.P)
        np.copyto(self.e_old, self.e)
        np.copyto(self.E_old, self.E)

    # Initialize variable with function at the spatial cell centers
    def initializeAtCenters(self, function):
        values = np.zeros(self.N)
        if function is not None:
            for i in range(self.N):
                values[i] = function(self.geo.r[i])
        return values

    # Initialize variable with function at the spatial cell edges
    def initializeAtEdges(self, function):
        values = np.zeros(self.N + 1)
        if function is not None:
            for i in range(self.N + 1):
                values[i] = function(self.geo.r_half[i])
        return values

    # Recompute radiation energy with updated internal energy
    def recomputeVelocity(self):
        self.u = self.momentum / self.rho

    # Recompute temperature with updated internal energy
    def recomputeInternalEnergy(self):
        self.e = self.mat_energy / self.rho  \
                 - 0.5 * (self.momentum / self.rho)**2

    # Recompute temperature with updated internal energy
    def recomputeT(self):
        self.T = (self.mat_energy / self.rho - \
                 0.5 * (self.momentum / self.rho)**2) \
                 / self.mat.C_v

    # Recompute pressure with updated density and internal energy
    def recomputeP(self):
        self.P = (self.mat.gamma - 1) * self.rho * 
                 self.mat_energy / self.rho  \
                 - 0.5 * (self.momentum / self.rho)**2
