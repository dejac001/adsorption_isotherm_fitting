"""Default solver for solving equations"""

import pyomo.environ as pyo
default_solver = pyo.SolverFactory('ipopt')
