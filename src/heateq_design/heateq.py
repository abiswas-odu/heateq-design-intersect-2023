from abc import ABC, abstractmethod
import numpy as np
import random
import math


class HeatEq(ABC):
    def __init__(self, lenx: float, maxt: float, alpha: float, dx: float,
                 dt: float, bc0: float, bc1: float, ic: str, outi: int):
        self.alpha = alpha
        self.dx = dx
        self.dt = dt
        self.bc0 = bc0
        self.bc1 = bc1
        self.ic = ic
        self.lenx = lenx
        self.maxt = maxt
        self.max_iter = 99999
        self.outi = outi

        self.Nx = int(self.lenx / self.dx) + 1
        self.Nt = int(self.maxt / self.dt)
        self.dx = self.lenx / (self.Nx - 1)

        # Init vectors
        self.curr = np.zeros(self.Nx)
        self.last = np.zeros(self.Nx)
        self.exact = np.zeros(self.Nx)
        self.change_history = np.zeros(self.Nx)
        self.error_history = np.zeros(self.Nx)

    def set_initial_condition(self):
        if self.ic.startswith("const("):  # const(val)
            cval = float(self.ic[self.ic.find("(") + 1: self.ic.find(")")])
            self.last[self.last == 0] = cval
        elif self.ic.startswith("step("):  # step(left,xmid,right)
            cvals = self.ic[self.ic.find("(") + 1: self.ic.find(")")].split(sep=',')
            left = float(cvals[0])
            xmid = float(cvals[1])
            right = float(cvals[2])
            x = 0
            for idx in range(self.Nx):
                if x < xmid:
                    self.last[idx] = left
                else:
                    self.last[idx] = right
                x = x + self.dx
        elif self.ic.startswith("ramp("):  # ramp(left,right)
            cvals = self.ic[self.ic.find("(") + 1: self.ic.find(")")].split(sep=',')
            left = float(cvals[0])
            right = float(cvals[1])
            dv = (right - left) / (self.Nx - 1)
            self.last[0] = left
            for idx in range(1, self.Nx):
                self.last[idx] = self.last[idx - 1] + dv
        elif self.ic.startswith("rand("):  # rand(seed,base,amp)
            cvals = self.ic[self.ic.find("(") + 1: self.ic.find(")")].split(sep=',')
            seed = int(cvals[0])
            base = float(cvals[1])
            amp = float(cvals[2])
            maxr = (1 << 31) - 1
            random.seed(seed)
            for idx in range(self.Nx):
                self.last[idx] = base + amp * (2 * random.random() / maxr - 1)
        elif self.ic.startswith("sin("):  # sin(PI*x)
            self.last[0] = self.dx
            for idx in range(1, self.Nx):
                self.last[idx] = math.pi * (self.last[idx - 1] + self.dx)
        elif self.ic.startswith("spikes("):  # spikes(Const,Amp,Loc,Amp,Loc,...)
            cvals = self.ic[self.ic.find("(") + 1: self.ic.find(")")].split(sep=',')
            self.last[self.last == 0] = float(cvals[0])
            for idx in range(1, len(cvals) - 1):
                amp = float(cvals[idx])
                spike_idx = int(cvals[idx + 1])
                if spike_idx < self.Nx:
                    self.last[spike_idx] = amp

    @abstractmethod
    def initialize(self):
        pass

    @abstractmethod
    def update_solution(self):
        pass

    def solve(self):
        self.initialize()

        # Iterate to max iterations or solution change is below threshold
        ti = 0
        while (ti * self.dt) < self.maxt:
            if not self.update_solution():
                print("Solution criteria violated. Make better choices\n")
                return

            # compute amount of change in solution
            diff = self.curr - self.last
            change = np.sum(diff*diff)

            # Handle possible termination by change threshold
            if self.maxt == self.max_iter or change < (-self.maxt * -self.maxt):
                print("Stopped after {0} iterations for threshold {1}\n".format(ti, change))
                return

            if self.outi and ti % self.outi == 0:
                print("Iteration {0}: last change l2={1}\n".format(ti, change))

            # Copy current solution to last
            self.last = self.curr
            ti = ti + 1
