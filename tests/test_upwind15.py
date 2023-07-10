from heateq_design.upwind15 import UpWind15
import numpy as np
from pytest import approx

def test_upwind():
    runame = 'upwind15_test_upwind'
    heat_solver = UpWind15(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 100, 10)
    heat_solver.solve(runame)
    cpp_result = np.array([0, 0.1978, 0.3835, 0.5473, 0.6832, 0.7892, 0.867, 0.9211, 0.9572, 0.9815, 1])
    assert heat_solver.curr == approx(cpp_result, rel=1e-3, abs=0)
