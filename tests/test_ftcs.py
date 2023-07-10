from heateq_design.ftcs import FTCS
import numpy as np
from pytest import approx

def test_ftcs():
    runame = 'ftcs_test_ftcs'
    heat_solver = FTCS(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 100, 10)
    heat_solver.solve(runame)
    cpp_result = [0, 0.1039, 0.2073, 0.3101, 0.4119, 0.5125, 0.6119, 0.7101, 0.8073, 0.9039, 1]
    assert heat_solver.curr == approx(cpp_result, rel=1e-3, abs=0)
