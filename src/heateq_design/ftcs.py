from .heateq import HeatEq


class FTCS(HeatEq):
    """
    This class represents the Forward-Time Central-Space (FTCS) scheme for solving the heat equation.

    Methods:
        initialize():
            Initializes the FTCS scheme by setting the initial conditions.

        update_solution():
            Updates the solution using the FTCS algorithm.

    Attributes:
        Inherits attributes from the base class HeatEq.

    """
    def initialize(self):
        """
        Initializes the FTCS scheme by setting the initial conditions.
        """
        self.set_initial_condition()

    def update_solution(self):
        """
        Updates the solution using the FTCS algorithm.

        Returns:
            bool: True if the update is successful and within stability limits, False otherwise.
        """
        r = self.alpha * self.dt / (self.dx * self.dx)

        # sanity check for stability
        if r > 0.5:
            return False

        # FTCS update algorithm
        for idx in range(1, self.Nx - 1):
            self.curr[idx] = r * self.last[idx + 1] + \
                             (1 - 2 * r) * self.last[idx] + \
                             r * self.last[idx - 1]

        # enforce boundary conditions
        self.curr[0] = self.bc0
        self.curr[self.Nx - 1] = self.bc1

        return True
