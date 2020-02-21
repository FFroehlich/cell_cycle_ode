import petab
import pypesto
import sys
import os
import numpy as np

import matplotlib.pyplot as plt

model_name = sys.argv[1]
cell_line = sys.argv[2]
agent = sys.argv[3]
regpar = np.float(sys.argv[4])

setting = f'{model_name}__{cell_line}__{agent}__lambda{regpar}'

n_threads = int(sys.argv[5])

basedir = os.path.dirname(os.path.realpath(__file__))
petab_problem = petab.Problem.from_yaml(
    os.path.join(basedir, 'petab',  f'{setting}.yaml')
)

importer = pypesto.PetabImporter(
    petab_problem,
    output_folder=os.path.join(basedir, 'models', model_name, model_name),
    model_name=model_name,
)

obj = importer.create_objective()
obj.n_threads = n_threads
obj.use_amici_petab_simulate = False
problem = importer.create_problem(obj)

optim_options = {
    'xtol': 1e-12,
    'gtol': 1e-4,
}

optimizer = pypesto.ScipyOptimizer(
    method='ls_trf',
    options=optim_options
)

optimize_options = pypesto.optimize.optimize.OptimizeOptions(
    startpoint_resample=True,
    allow_failed_starts=True,
)

result = pypesto.minimize(
    problem=problem,
    optimizer=optimizer,
    n_starts=10,
    options=optimize_options
)

resultdir = os.path.join(basedir, 'results')
os.makedirs(resultdir, exist_ok=True)
pypesto.visualize.parameters(result, start_indices=range(5))
plt.savefig(
    os.path.join(resultdir, f'parameters__{setting}.pdf')
)

pypesto.visualize.waterfall(result, start_indices=range(5))
plt.savefig(
    os.path.join(resultdir, f'waterfall__{setting}.pdf')
)

result.optimize_result.as_dataframe().to_csv(
    os.path.join(resultdir, f'multistarts__{setting}.csv')
)
