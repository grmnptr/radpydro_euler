from inputparameters import InputParameters
from radpydro import RadPydro
import numpy as np
import matplotlib.pyplot as plt

input = InputParameters()
input.geometry = 'slab'
input.N = 5
input.r_half = np.linspace(0, 1, num=input.N + 1) # cm
input.C_v = 1.66 # jerks / (cm3 eV)
input.gamma = 1.5 # cm3 / g
input.kappa = [1, 1, 1, 1] # g/cm2
input.kappa_s = 1 # g / cm2
input.a = 0.01372 # [jerks / (cm3 kev4)]
input.c = 299.792 # [cm / sh]

# Initial conditions
input.rho = lambda r: 1.0 #g/cm3
input.T = lambda r: 1.0 # keV
input.u = lambda r: 0.0 # cm / sh
input.E = lambda r: input.a * input.T(0)**4

# Boundary conditions
input.hydro_L = 'u'
input.hydro_L_val = 0.01
input.hydro_R = 'u'
input.hydro_R_val = 0.01
input.rad_L = 'source'
input.rad_L_val = input.a * input.T(0)**4
input.rad_R = 'source'
input.rad_R_val = 0

# Iteration controls
input.CoFactor = 1
input.relEFactor = 0.2
input.maxTimeStep = 0.01
input.T_final = 1

# Slope control
input.omega = 0

rp = RadPydro(input)
rp.computeTimeStep()