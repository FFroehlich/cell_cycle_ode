import os
import yaml

import pandas as pd
import numpy as np

from cc_model import model

petab_dir = 'petab'

df = pd.read_csv(os.path.join('data', 'cell_cycle_data.csv'))

observables = pd.DataFrame({
    'observableId': [obs.name for obs in model.observables],
    'observableName': [obs.name.replace('_obs', '')
                       for obs in model.observables],
    'observableFormula': [f'__obs{iobs}'
                          for iobs, obs in enumerate(model.observables)],
    'observableParameters': ['' for _ in model.observables],
    'noiseFormula': [obs.name for obs in model.observables],
})

observable_file = f'observables.tsv'
observables.to_csv(os.path.join(petab_dir, observable_file), sep='\t')

for line in df.cell_line.unique():

    df_line = df.query(f'cell_line == "{line}"').copy()

    measurement_cols = ['G1', 'S', 'G2_plus_M', 'D']
    df_line.dropna(axis=0, how='all', subset=measurement_cols, inplace=True)

    measurements = pd.melt(df_line,
                           id_vars=['timepoint', 'agent',
                                    'concentration'],
                           value_vars=measurement_cols,
                           value_name='measurement',
                           var_name='observableId')

    measurements.rename(columns={
        'timepoint': 'time',
    }, inplace=True)

    measurements.observableId = measurements.observableId.apply(
        lambda x: f'{x}_obs'
    )

    measurements.dropna(axis=0,  subset=['measurement'], inplace=True)

    measurements['simulationConditionId'] = \
        measurements[['agent', 'concentration']].apply(
            lambda x: f'{x[0]}__{np.log10(x[1]):.2}'
            .replace('__nan', '').replace('.', '_'),
            axis=1
        )

    measurements['observableParameters'] = ''
    measurements['noiseParameters'] = ''

    measurement_file = f'{line}_measurements.tsv'
    measurements[['observableId', 'time', 'measurement', 'noiseParameters',
                  'simulationConditionId', 'observableParameters']].to_csv(
        os.path.join(petab_dir, measurement_file), sep='\t', index=False
    )

    conditions = pd.DataFrame([
        {
            'conditionId': cond,
            'rphi': f'{cond}_rphi' if cond != 'DMSO' else 1.0,
            'rG1S': f'{cond}_rG1S' if cond != 'DMSO' else 1.0,
            'rSG2': f'{cond}_rSG2' if cond != 'DMSO' else 1.0,
            'rG2MG1': f'{cond}_rG2MG1' if cond != 'DMSO' else 1.0,
        }
        for cond in measurements['simulationConditionId'].unique()
    ])
    condition_file = f'{line}_conditions.tsv'
    conditions.to_csv(os.path.join(petab_dir, condition_file),
                      sep='\t', index=False)

    baseline_pars = ['G1_0', 'S_0', 'G2_0', 'D_0', 'kphi', 'kG1S', 'kSG2',
                     'kG2MG1']
    condition_pars = np.unique(
        conditions.loc[conditions.conditionId != 'DMSO',
                       ['rphi', 'rSG2', 'rG1S', 'rG2MG1']].values
    )

    parameters = pd.concat([
        pd.DataFrame(
            {
                'parameterId': baseline_pars,
                'parameterScale': ['log10' for _ in baseline_pars],
                'lowerBound': [1e-5 for _ in baseline_pars],
                'upperBound': [10 for _ in baseline_pars],
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
                #'prior': ['logLaplace' for _ in condition_pars]
            }
        ),
    ])
    parameter_file = f'{line}_parameters.tsv'
    parameters.to_csv(os.path.join(petab_dir, parameter_file),
                      sep='\t', index=False)

    petab_info = {
        'format_version': 0,
        'parameter_file': parameter_file,
        'problems': [
            {
                'sbml_files': ['cell_cycle_no_M.sbml'],
                'measurement_files': [measurement_file],
                'condition_files': [condition_file],
                'observable_files': [observable_file],
            }
        ]
    }

    with open(os.path.join(petab_dir, f'{line}.yaml'), 'w') as file:
        yaml.dump(petab_info, file)
