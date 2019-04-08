import numpy as np
from sys import exit

def isScalar(value):
    return isinstance(value, float) or isinstance(value, int)

class InputParameters:
    def __init__(self):
        # Geometry specficiations
        self.geometry = 'slab'
        self.r_half = None

        # Whether or not to enable the hydro/radiation run
        self.enable_hydro = True
        self.enable_radiation = True

        # Material properties
        self.C_v = None
        self.gamma = None
        self.kappa = None
        self.kappa_s = None

        # Constants
        self.a = None
        self.c = None

        # Initial conditions
        self.E = None
        self.rho = None
        self.T = None
        self.u = None

        # Left ghost cell
        self.rho_bL = None
        self.e_bL = None
        self.u_bL = None

        # Right ghost cell
        self.rho_bR = None
        self.e_bR = None
        self.u_bR = None

        # Iteration parameters
        self.CoFactor = None
        self.relEFactor = None
        self.maxTimeStep = None
        self.T_final = None

        # Slope withing factor
        self.omega = None
