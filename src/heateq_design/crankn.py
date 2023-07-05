from .heateq import HeatEq
import numpy as np


class CrankN(HeatEq):
    """
    This class represents the Crank-Nicolson scheme for solving the heat equation.

    Args:
        lenx (float): Length of the domain in the x-direction.
        maxt (float): Maximum time for the simulation.
        alpha (float): Thermal diffusivity.
        dx (float): Spatial step size.
        dt (float): Time step size.
        bc0 (float): Boundary condition at x = 0.
        bc1 (float): Boundary condition at x = lenx.
        ic (str): Initial condition type.
        outi (int): Output interval.

    Methods:
        __init__(lenx, maxt, alpha, dx, dt, bc0, bc1, ic, outi):
            Initializes the Crank-Nicolson scheme.

        r83_np_fa():
            Factors the tridiagonal matrix.

        initialize():
            Initializes the Crank-Nicolson scheme by setting the initial conditions.

        update_solution():
            Updates the solution using the Crank-Nicolson algorithm.

    Attributes:
        Inherits attributes from the base class HeatEq.

    """
    def __init__(self, lenx: float, maxt: float, alpha: float, dx: float, dt: float, bc0: float, bc1: float, ic: str,
                 outi: int):
        """
        Initializes the Crank-Nicolson scheme.

        Args:
            lenx (float): Length of the domain in the x-direction.
            maxt (float): Maximum time for the simulation.
            alpha (float): Thermal diffusivity.
            dx (float): Spatial step size.
            dt (float): Time step size.
            bc0 (float): Boundary condition at x = 0.
            bc1 (float): Boundary condition at x = lenx.
            ic (str): Initial condition type.
            outi (int): Output interval.
        """
        super().__init__(lenx, maxt, alpha, dx, dt, bc0, bc1, ic, outi)
        w = self.alpha * self.dt / self.dx / self.dx

        # Build a tri-diagonal matrix
        self.cn_Amat = np.zeros(3 * self.Nx)

        self.cn_Amat[0 + 0 * 3] = 0.0
        self.cn_Amat[1 + 0 * 3] = 1.0
        self.cn_Amat[0 + 1 * 3] = 0.0

        for idx in range(1, self.Nx - 1):
            self.cn_Amat[2 + (idx - 1) * 3] = - w
            self.cn_Amat[1 + idx * 3] = 1.0 + 2.0 * w
            self.cn_Amat[0 + (idx + 1) * 3] = - w

        self.cn_Amat[2 + (self.Nx - 2) * 3] = 0.0
        self.cn_Amat[1 + (self.Nx - 1) * 3] = 1.0
        self.cn_Amat[2 + (self.Nx - 1) * 3] = 0.0

        # Factor the matrix.
        self.r83_np_fa()

    def r83_np_fa(self):
        """
        Factors the tridiagonal matrix.
        """
        for idx in range(1, self.Nx):
            assert (self.cn_Amat[1 + (idx - 1) * 3] != 0.0)

            # Store the multiplier in L.
            self.cn_Amat[2 + (idx - 1) * 3] = self.cn_Amat[2 + (idx - 1) * 3] / \
                                              self.cn_Amat[1 + (idx - 1) * 3]

            # Modify the diagonal entry in the next column.
            self.cn_Amat[1 + idx * 3] = self.cn_Amat[1 + idx * 3] - \
                                        self.cn_Amat[2 + (idx - 1) * 3] * \
                                        self.cn_Amat[0 + idx * 3]
        assert (self.cn_Amat[1 + (self.Nx - 1) * 3] != 0.0)

    def initialize(self):
        """
        Factors the tridiagonal matrix.
        """
        self.set_initial_condition()

    def update_solution(self):
        """
        Updates the solution using the Crank-Nicolson algorithm.

        Returns:
            bool: True if the update is successful, False otherwise.
        """
        # r83_np_sl
        for idx in range(self.Nx):
            self.curr[idx] = self.last[idx]

        # Solve L * Y = B.
        for idx in range(1, self.Nx):
            self.curr[idx] = self.curr[idx] - \
                             self.cn_Amat[2+(idx-1)*3] * \
                             self.curr[idx-1]

        # Solve U * X = Y.
        for idx in range(self.Nx, 0, -1):
            self.curr[idx-1] = self.curr[idx-1] / self.cn_Amat[1+(idx-1)*3]
            if 1 < idx:
                self.curr[idx-2] = self.curr[idx-2] - \
                                   self.cn_Amat[0+(idx-1)*3] * \
                                   self.curr[idx-1]

        self.curr[0] = self.bc0
        self.curr[self.Nx - 1] = self.bc1
        return True
