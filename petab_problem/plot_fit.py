import sys
import os

import petab
from petab.visualize import plot_data_and_simulation
import pypesto
from plotnine import *
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np


model_name = sys.argv[1]
cell_line = sys.argv[2]
agent = sys.argv[3]
regpar = np.float(sys.argv[4])

setting = f'{model_name}__{cell_line}__{agent}__lambda{regpar}'

basedir = os.path.dirname(os.path.realpath(__file__))
resultdir = os.path.join(basedir, 'results')
results = pd.read_csv(os.path.join(resultdir,
                                   f'multistarts__{setting}.csv'))

petab_problem = petab.Problem.from_yaml(
    os.path.join(basedir, 'petab',  f'{setting}.yaml')
)

value_name = 'value'
x = np.fromstring(results.loc[0, 'x'][1:-2], sep=' ')

importer = pypesto.PetabImporter(
    petab_problem,
    output_folder=os.path.join(basedir, 'models', model_name, model_name),
    model_name=model_name,
)

obj = importer.create_objective()
ret = obj(x, return_dict=True)

simu_meas_df = importer.rdatas_to_measurement_df(ret['rdatas'])

simu_meas_df = simu_meas_df.rename(columns={"measurement": "simulation"})


vis_spec = pd.DataFrame({
    'plotId':
        ['control', 'control', 'control', 'control',
         'perturbation1', 'perturbation2', 'perturbation3', 'perturbation4'],
    'plotName':
        ['DMSO', 'DMSO', 'DMSO', 'DMSO',
         agent, agent, agent, agent],
    'plotTypeSimulation':
        ['LinePlot', 'LinePlot', 'LinePlot', 'LinePlot',
         'LinePlot', 'LinePlot', 'LinePlot', 'LinePlot'],
    'plotTypeData':
        ['MeanAndSD', 'MeanAndSD', 'MeanAndSD', 'MeanAndSD',
         'MeanAndSD', 'MeanAndSD', 'MeanAndSD', 'MeanAndSD'],
    'datasetId':
        ['control', 'control', 'control', 'control',
         'perturbation', 'perturbation', 'perturbation', 'perturbation'],
    'xValues':
        ['time', 'time', 'time', 'time',
         'concentration', 'concentration', 'concentration', 'concentration'],
    'xScale':
        ['lin', 'lin', 'lin', 'lin',
         'log10', 'log10', 'log10', 'log10'],
    'xLabel':
        ['time', 'time', 'time', 'time',
         'concentration', 'concentration', 'concentration', 'concentration'],
    'yValues':
        ['D_obs', 'S_obs', 'G1_obs', 'G2_obs',
         'D_obs', 'S_obs', 'G1_obs', 'G2_obs'],
    'yLabel':
        ['cell count', 'cell count', 'cell count', 'cell count',
         'cell count', 'cell count', 'cell count', 'cell count'],
    'yScale': [
        'log10', 'log10', 'log10', 'log10',
        'log10', 'log10', 'log10', 'log10'],
    'legendEntry':
        ['D', 'S', 'G1', 'G2_plus_M',
         'D', 'S', 'G1', 'G2_plus_M'],
    'xOffset':
        [-3, -1, 1, 3,
         0.0, 0.0, 0.0, 0.0],
})

plot_data_and_simulation(
    exp_data=petab_problem.measurement_df,
    exp_conditions=petab_problem.condition_df,
    sim_data=simu_meas_df,
    vis_spec=vis_spec,
)
plt.savefig(os.path.join(resultdir, f'simulation_and_data__{setting}.pdf'))
