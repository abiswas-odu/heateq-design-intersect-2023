from heateq_design.crankn import CrankN
from pytest import approx
import numpy as np

def test_1():
    runame = 'crankn_test_1'
    heat_solver = CrankN(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 100, 10)
    heat_solver.solve(runame)
    cpp_result = np.array([0, 0.104, 0.2076, 0.3104, 0.4122, 0.5129, 0.6122, 0.7104, 0.8076, 0.904, 1])
    assert heat_solver.curr == approx(cpp_result, rel=1e-3, abs=0)
