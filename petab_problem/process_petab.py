import os
import yaml

import pandas as pd
import numpy as np

import sys
import importlib
import petab


def get_measurements_df(df, model, cols):
    df.dropna(axis=0, how='all', subset=cols, inplace=True)

    measurements = pd.melt(df,
                           id_vars=['timepoint',
                                    'agent',
                                    'concentration'],
                           value_vars=cols,
                           value_name='measurement',
                           var_name='observableId')

    measurements.rename(columns={
        'timepoint': 'time',
    }, inplace=True)

    measurements.observableId = measurements.observableId.apply(
        lambda x: f'{x}_obs'
    )

    measurements.dropna(axis=0, subset=['measurement'], inplace=True)

    measurements['simulationConditionId'] = \
        measurements[['agent', 'concentration']].apply(
            lambda x: f'{x[0]}__{np.log10(x[1]):.2}'.replace('.', '_')
            if x[0] != 'DMSO' else 'DMSO',
            axis=1
        )

    measurements['datasetId'] = measurements['simulationConditionId'].apply(
        lambda x: 'control' if x == 'DMSO' else 'perturbation'
    )

    measurements['observableParameters'] = ''
    measurements['preequilibrationConditionId'] = ''
    measurements['noiseParameters'] = ''

    for simCond in measurements.simulationConditionId.unique():
        subset = measurements.query(f'simulationConditionId == "{simCond}"')
        for time in subset.time.unique():
            for obs in subset.query(f'time == {time}').observableId.unique():
                sel = (measurements.time == time) \
                      & (measurements.observableId == obs) \
                      & (measurements.simulationConditionId == simCond)
                measurements.loc[sel, 'noiseParameters'] = \
                    np.max([measurements.loc[sel, 'measurement'].std(), 1.0])

        # add regularization term for each non-DMSO condition
        for expr in model.expressions:
            if not expr.name.endswith('_obs'):
                continue
            measurements = measurements.append({
                'time': 0.0,
                'agent': '',
                'concentration': 0.0,
                'observableId': expr.name,
                'measurement': 0.0,
                'simulationConditionId': simCond,
                'preequilibrationConditionId': '',
                'observableParameters': '',
                'noiseParameters': '',
                'datasetId': 'regularization'
            }, ignore_index=True)

    measurements.drop(columns=['agent'], inplace=True)
    measurements.sort_values('concentration', inplace=True)
    return measurements


def get_conditions_df(measurements):
    conditions = pd.DataFrame([
        {
            'conditionId': cond,
            'rphi': f'{cond}_rphi' if cond != 'DMSO' else 1.0,
            'rG1S': f'{cond}_rG1S' if cond != 'DMSO' else 1.0,
            'rSG2': f'{cond}_rSG2' if cond != 'DMSO' else 1.0,
            'rG2MG1': f'{cond}_rG2MG1' if cond != 'DMSO' else 1.0,
            'concentration': measurements.loc[
                (measurements['simulationConditionId'] == cond)
                & (measurements['observableId'] == 'D_obs'), # avoid reg
                'concentration'].values[0]
        }
        for cond in measurements['simulationConditionId'].unique()
    ])
    conditions.sort_values('concentration', inplace=True)
    return conditions


def get_parameters_df(conditions, baseline_pars, condition_ratios, regpar):
    condition_pars = np.unique(
        conditions.loc[conditions.conditionId != 'DMSO',
                       condition_ratios].values
    )

    parameters = pd.concat([
        pd.DataFrame(
            {
                'parameterId': baseline_pars,
                'parameterScale': ['log10' for _ in baseline_pars],
                'lowerBound': [
                    1 if par.endswith('_0') else 1e-5
                    for par in baseline_pars
                ],
                'upperBound': [
                    1e4 if par.endswith('_0') else 10
                    for par in baseline_pars
                ],
                'estimate': [True for _ in baseline_pars],
                'nominalValue': ['' for _ in baseline_pars]
            }
        ),
        pd.DataFrame(
            {
                'parameterId': list(condition_pars),
                'parameterScale': ['log10' for _ in condition_pars],
                'lowerBound': [1e-3 for _ in condition_pars],
                'upperBound': [1e3 for _ in condition_pars],
                'estimate': [True for _ in condition_pars],
                'nominalValue': ['' for _ in condition_pars]
                # 'prior': ['logLaplace' for _ in condition_pars]
            }
        ),
        pd.DataFrame(
            {
                'parameterId': ['lambda_reg'],
                'parameterScale': ['log10'],
                'lowerBound': [1e-3],
                'upperBound': [1e3],
                'estimate': [False],
                'nominalValue': [regpar]
            }
        )
    ])
    return parameters


def process_to_petab(df, model, regpar, prefix):
    measurements = get_measurements_df(
        df,
        model,
        ['G1', 'S', 'G2_plus_M', 'D']
    )
    conditions = get_conditions_df(measurements)

    measurements.drop(columns=['concentration'], inplace=True)

    condition_pars = ['rphi', 'rSG2', 'rG1S', 'rG2MG1']
    baseline_pars = [par.name for par in model.parameters
                     if par.name not in condition_pars]

    parameters = get_parameters_df(conditions,
                                   baseline_pars,
                                   condition_pars,
                                   regpar)

    measurement_file = f'{prefix}_measurements.tsv'
    condition_file = f'{prefix}_conditions.tsv'
    parameter_file = f'{prefix}_parameters.tsv'

    file_data = {
        measurement_file: measurements,
        condition_file: conditions,
        parameter_file: parameters,
    }
    for file, data in file_data.items():
        data.to_csv(os.path.join(petab_dir, file),
                    sep='\t', index=False)

    return measurement_file, parameter_file, condition_file


def export_petab_yaml(parameter_file, measurement_files, condition_files,
                      observable_file, model_name, name,):
    petab_info = {
        'format_version': petab.format_version.__format_version__,
        'parameter_file': parameter_file,
        'problems': [
            {
                'sbml_files': [f'{model_name}.sbml'],
                'measurement_files': measurement_files,
                'condition_files': condition_files,
                'observable_files': [observable_file],
            }
        ]
    }

    with open(os.path.join(petab_dir, f'{name}.yaml'), 'w') as file:
        yaml.dump(petab_info, file)


model_name = sys.argv[1]
cell_line = sys.argv[2]
agent = sys.argv[3]
regpar = np.float(sys.argv[4])

setting = f'{model_name}__{cell_line}__{agent}__lambda{regpar}'

basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(basedir, 'pysb'))
model = importlib.import_module(f'model_{model_name}').model

petab_dir = os.path.join(basedir, 'petab')
os.makedirs(petab_dir, exist_ok=True)

df = pd.read_csv(os.path.join('data', 'cell_cycle_data.csv'))

# fix misannotation
df.loc[df.agent == 'Storausporin', 'agent'] = 'Staurosporine'

if not sum(df.cell_line == cell_line):
    raise ValueError(f'Unknownt cell line {cell_line}')

if not sum(df.agent == agent):
    raise ValueError(f'Unknownt agent {agent}')

observables = pd.DataFrame({
    'observableId': [obs.name for obs in model.observables],
    'observableName': [obs.name.replace('_obs', '')
                       for obs in model.observables],
    'observableFormula': [f'__obs{iobs}'
                          for iobs, obs in enumerate(model.observables)],
    'observableParameters': ['' for _ in model.observables],
    'noiseFormula': [f'noiseParameter1_{obs.name}'
                     for obs in model.observables],
    #'observableTransformation': ['log10' for obs in model.observables],
})

for expr in model.expressions:
    if not expr.name.endswith('obs'):
        continue

    observables = observables.append({
        'observableId': expr.name,
        'observableName': f'ridge regularization {expr.name.split("_")[0]}',
        'observableFormula': expr.name,
        'observableParameters': '',
        'noiseFormula': '1/sqrt(lambda_reg)',
    }, ignore_index=True)

observable_file = f'{setting}_observables.tsv'
observables.to_csv(os.path.join(petab_dir, observable_file), sep='\t',
                   index=False)

df_line = df.query(f'cell_line == "{cell_line}"')

df_line_dmso = df_line.query(f'agent == "DMSO"').copy()
df_line_agent = df_line.query(f'agent == "{agent}"').copy()

m_file, p_file, c_file = process_to_petab(
    pd.concat([df_line_agent, df_line_dmso]),
    model,
    regpar,
    setting
)

export_petab_yaml(p_file,
                  [m_file],
                  [c_file],
                  observable_file,
                  model_name,
                  setting)


