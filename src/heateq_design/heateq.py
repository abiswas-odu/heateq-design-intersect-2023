import os.path
from abc import ABC, abstractmethod
import numpy as np
import random
import math
import shutil


def write_array(file_name, var_name, dx, a):
    with open(file_name, 'w') as out_f:
        out_f.write('# {0}\n'.format(var_name))
        for i in range(0, len(a)):
            out_f.write('{0} {1}\n'.format(i*dx, a[i]))


class HeatEq(ABC):
    """
    Abstract base class representing the heat equation.

    Methods:
        set_initial_condition():
            Sets the initial condition based on the specified string.

        initialize():
            Abstract method to be implemented by subclasses. Initializes the heat equation.

        update_solution():
            Abstract method to be implemented by subclasses. Updates the solution.

        solve():
            Solves the heat equation by iterating until the maximum number of iterations or a change threshold is reached.

    Attributes:
        lenx (float): Length of the domain.
        maxt (float): Maximum time.
        alpha (float): Thermal diffusivity.
        dx (float): Spatial step size.
        dt (float): Time step size.
        bc0 (float): Boundary condition at x = 0.
        bc1 (float): Boundary condition at x = lenx.
        ic (str): Initial condition string.
        outi (int): Output interval.
        max_iter (int): Maximum number of iterations.
        Nx (int): Number of spatial grid points.
        Nt (int): Number of time steps.
        curr (np.ndarray): Current solution vector.
        last (np.ndarray): Solution vector from the previous time step.
        exact (np.ndarray): Exact solution vector.
        change_history (np.ndarray): History of solution changes.
        error_history (np.ndarray): History of solution errors.

    """

    def __init__(self, lenx: float, maxt: float, alpha: float, dx: float,
                 dt: float, bc0: float, bc1: float, ic: str, outi: int):
        """
        Initializes the HeatEq class with the specified parameters.

        Args:
            lenx (float): Length of the domain.
            maxt (float): Maximum time.
            alpha (float): Thermal diffusivity.
            dx (float): Spatial step size.
            dt (float): Time step size.
            bc0 (float): Boundary condition at x = 0.
            bc1 (float): Boundary condition at x = lenx.
            ic (str): Initial condition string.
            outi (int): Output interval.

        """
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
        """
        Sets the initial condition based on the specified string.
        """
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
        """
        Abstract method to be implemented by subclasses.
        Initializes the heat equation.
        """
        pass

    @abstractmethod
    def update_solution(self):
        """
        Abstract method to be implemented by subclasses.
        Updates the solution.
        """
        pass

    def solve(self, output_name):
        """
        Solves the heat equation by iterating until the maximum number of iterations or a change threshold is reached.
        """
        if os.path.isdir(output_name):
            shutil.rmtree(output_name)
        os.makedirs(output_name)

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
                break

            if self.outi and ti % self.outi == 0:
                print("Iteration {0}: last change l2={1}\n".format(ti, change))

            # Copy current solution to last
            self.last = self.curr
            ti = ti + 1
        write_array(os.path.join(output_name, output_name + '_soln_final.curve'),
                         'Temperature', self.dx, self.curr)

