from src.heateq_design.ftcs import *

def test_ftcs_1():
    runame = 'ftcs_test_1'
    heat_solver = FTCS(1.0, 2.0, 0.2, 0.1, 0.004, 0, 1, 'const(1)', 1, 1)
    heat_solver.solve(runame)