from .heateq import HeatEq


class UpWind15(HeatEq):
    """
    This class represents the Upwind 1.5 scheme for solving the heat equation.

    Methods:
        initialize():
            Initializes the Upwind 1.5 scheme by setting the initial conditions.

        update_solution():
            Updates the solution using the Upwind 1.5 algorithm.

    Attributes:
        Inherits attributes from the base class HeatEq.

    """
    def initialize(self):
        """
        Initializes the Upwind 1.5 scheme by setting the initial conditions.
        """
        self.set_initial_condition()

    def update_solution(self):
        """
        Updates the solution using the Upwind 1.5 algorithm.

        Returns:
            bool: True if the update is successful, False otherwise.
        """
        f2 = 1.0 / 24
        f1 = 1.0 / 6
        f0 = 1.0 / 4
        k = self.alpha * self.alpha * self.dt / (self.dx * self.dx)
        k2 = k * k

        self.curr[0] = self.bc0
        self.curr[1] = self.last[1] + k * (self.last[0] - 2 * self.last[1] + self.last[2])
        self.curr[self.Nx - 2] = self.last[self.Nx - 2] + k * \
                                 (self.last[self.Nx - 3] -
                                  2 * self.last[self.Nx - 2] +
                                  self.last[self.Nx - 1])
        self.curr[self.Nx - 1] = self.bc1
        for idx in range(2, self.Nx - 2):
            self.curr[idx] = f2 * (12 * k2 - 2 * k) * self.last[idx - 2] \
                             + f2 * (12 * k2 - 2 * k) * self.last[idx + 2] \
                             - f1 * (12 * k2 - 8 * k) * self.last[idx - 1] \
                             - f1 * (12 * k2 - 8 * k) * self.last[idx + 1] \
                             + f0 * (12 * k2 - 10 * k + 4) * self.last[idx]

        return True
