from src.heateq_design.upwind15 import *
import numpy as np

def test_1():
    runame = 'upwind15_test_1'
    heat_solver = UpWind15(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 100, 10)
    heat_solver.solve(runame)
    cpp_result = [0, 0.1978, 0.3835, 0.5473, 0.6832, 0.7892, 0.867, 0.9211, 0.9572, 0.9815, 1]
    np.testing.assert_allclose(heat_solver.curr, cpp_result, rtol=1e-3, atol=0)
