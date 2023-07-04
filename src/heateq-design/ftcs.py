from heateq import HeatEq


class FTCS(HeatEq):
    def initialize(self):
        self.set_initial_condition()

    def update_solution(self):
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
