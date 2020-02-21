import pandas as pd
import numpy as np
import petab

import sys
import os

import matplotlib.pyplot as plt
import seaborn as sns


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
df_pars = pd.DataFrame(dict(zip(
    petab_problem.x_ids, np.fromstring(results.loc[0, 'x'][1:-2], sep=' ')
)), index=[value_name])

cond_pars = df_pars[[
    col
    for col in df_pars.columns
    if col.startswith(agent)
]].transpose()

cond_pars[value_name] = cond_pars[value_name].apply(lambda x: np.power(10, x))

cond_pars['cell_line'] = cell_line
cond_pars['agent'] = agent
cond_pars['index'] = cond_pars.index

cond_pars['parameter'] = cond_pars['index'].apply(
    lambda x: f'$k_{{{x.split("_")[-1][1:]}}}$'
)
conc_name = 'conc'
cond_pars[conc_name] = cond_pars['index'].apply(
    lambda x: np.power(10, np.float('.'.join(x.split('_')[-3:-1])))
)

sns.lineplot(x=conc_name,
             y=value_name,
             hue='parameter',
             markers='o',
             data=cond_pars,)
plt.title(f'{cell_line}')
plt.xlabel(f'concentration {agent} [$\mu M$]')
plt.ylabel('parameter deviation')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-3, 1e3])
plt.savefig(os.path.join(resultdir, f'rate_deviation__{setting}.pdf'))
