from src.heateq_design.ftcs import *
import numpy as np

def test_1():
    runame = 'ftcs_test_1'
    heat_solver = FTCS(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 100, 10)
    heat_solver.solve(runame)
    cpp_result = [0, 0.1039, 0.2073, 0.3101, 0.4119, 0.5125, 0.6119, 0.7101, 0.8073, 0.9039, 1]
    np.testing.assert_allclose(heat_solver.curr, cpp_result, rtol=1e-3, atol=0)