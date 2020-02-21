import sys
import os
import amici.petab_import
import petab.problem
import importlib
import numpy as np

import pysb.export


model_name = sys.argv[1]
cell_line = sys.argv[2]
agent = sys.argv[3]
regpar = np.float(sys.argv[4])

setting = f'{model_name}__{cell_line}__{agent}__lambda{regpar}'

basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(basedir, 'pysb'))
model = importlib.import_module(f'model_{model_name}').model

sbml_model = os.path.join(basedir, 'petab', f'{model_name}.sbml')

with open(sbml_model, 'w') as f:
    f.write(pysb.export.export(model, 'sbml'))

modeldir = os.path.join(basedir, 'models')
os.makedirs(modeldir, exist_ok=True)

amici.petab_import.import_petab_problem(
    petab_problem=petab.problem.Problem.from_yaml(
        os.path.join(basedir, 'petab',
                     f'{setting}.yaml')
    ),
    model_output_dir=os.path.join(modeldir, model_name),
    model_name=model_name,
    force_compile=True,
)
