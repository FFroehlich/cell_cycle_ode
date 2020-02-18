import sys
import os
import amici
import amici.petab_import
import petab
import importlib

import pysb.export

model_name = sys.argv[1]

model = importlib.import_module(model_name).model

sbml_model = os.path.join('petab', f'{model.name}.sbml')

with open(sbml_model, 'w') as f:
    f.write(pysb.export.export(model, 'sbml'))

amici.petab_import.import_petab_problem(
    petab_problem=,
    model_output_dir=model.name,
    model_name=model.name,
    force_compile=True,
)