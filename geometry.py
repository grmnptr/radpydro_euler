import numpy as np

class Geometry:
    def __init__(self, rp):
        self.rp = rp
        self.input = rp.input

        # Radius at cell edges (user defined)
        self.r_half = np.copy(self.input.r_half)

        # Number of cells
        self.N = self.r_half.size - 1

        # Areas (defined at edges)
        self.A = np.zeros(self.N + 1)

        # Volumes (defined on spatial cells)
        self.V = np.zeros(self.N)

        # Radius at cell centers
        self.r = np.zeros(self.N)

        # Cells widths
        self.dr = np.zeros(self.N)

        # Initialize A, V, r
        self.computeGeometry()
        self.stepGeometry()

    # Compute geometry 
    def recomputeGeometry(self):

        A = self.A
        V = self.V
        r_half = self.r_half
        r = self.r
        dr = self.dr

        # Update areas
        self.recomputeAreas(A, r_half)

        # Update volumes
        self.recomputeVolumes(V, r_half)
        
        # Update cell centered radii and cell widths
        for i in range(self.N):
            r[i] = (r_half[i] + r_half[i + 1]) / 2
            dr[i] = (r_half[i + 1] - r_half[i])

class SlabGeometry(Geometry):
    def recomputeAreas(self, A, r_half):
        A.fill(1.0)

    def recomputeVolumes(self, V, r_half):
        for i in range(self.N):
            V[i] = r_half[i + 1] - r_half[i]

class CylindricalGeometry(Geometry):
    def recomputeAreas(self, A, r_half):
        for i in range(self.N + 1):
            A[i] = 2 * np.pi * r_half[i]

    def recomputeVolumes(self, V, r_half):
        for i in range(self.N):
            V[i] = np.pi * (r_half[i + 1]**2 - r_half[i]**2)

class SphericalGeometry(Geometry):
    def recomputeAreas(self, A, r_half):
        for i in range(self.N + 1):
            A[i] = 4 * np.pi * r_half[i]**2

    def recomputeVolumes(self, V, r_half):
        for i in range(self.N):
            V[i] = 4 * np.pi * (r_half[i + 1]**3 - r_half[i]**3) / 3
