from model_cell_cycle_no_M import model

from sympy import log
from pysb import Expression, Model

model = Model('cell_cycle_no_M_ridge', base=model)

for par in model.parameters:
    if par.name.startswith('r'):
        Expression(f'{par.name}_obs', log(par))

