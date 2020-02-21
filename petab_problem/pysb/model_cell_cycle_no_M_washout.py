from model_cell_cycle_no_M import model
from pysb import Model, Parameter, Rule

model = Model('cell_cycle_no_M_washout', base=model)

Rule('washout', model.monomers.dead_cell() >> None, Parameter('kwash', 0.1))
