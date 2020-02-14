# ODE model for cell cycle

import numpy as np
import pandas as pd
import math
import itertools
from scipy.integrate import odeint
from scipy.optimize import minimize
from random import  seed, randrange, uniform
import statistics
import matplotlib
import  seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
#cell cycle phases for which the ratets need to be estimated
#phases=['G1', 'S', 'G2_plus_M', 'D']

def ode_soln(x0, t, logk):
    '''
    Solution to the cell cycle ODE's, initial condition: x(0) = x0
    '''
    def cell_cycle(x, t, logk):
        '''
        ODE's for the cell cycle transition model
        '''

        G1, S, G2, M, D = x                            
        # number of cells in G1, S, G2, M phases and dead cells
        kG1S, kSG2, kG2M, kMG1, kphi =\
        pow(10,logk[0]), pow(10,logk[1]), pow(10,logk[2]), pow(10,logk[3]), pow(10, logk[4])  
        # tansition rates: G1->S, S->G2, G2->M, (G1, S, G2, M)->D
        dxdt =  [-(kG1S + kphi)*G1 + 2.0*kMG1*M,
                -(kSG2 + kphi)*S + kG1S*G1,
                -(kG2M + kphi)*G2 + kSG2*S,
                -(kMG1 + kphi)*M + kG2M*G2,
                (G1 + S + G2 + M)*kphi]
        return dxdt
    return  odeint(cell_cycle, x0, t, args=(logk,), 
                   rtol=1e-4,  atol=1e-4, hmin=1e-4)


def ode_soln_noM(x0, t, logk):
    '''
    Solution to the cell cycle ODE's, initial condition: x(0) = x0
    '''
    def cell_cycle_noM(x, t, logk):
        '''
        ODE's for the cell cycle transition model
        '''

        G1, S, G2_plus_M, D = x                            
        # number of cells in G1, S, G2, M phases and dead cells
        kG1S, kSG2M, kG2MG1, kphi =\
        pow(10,logk[0]), pow(10,logk[1]), pow(10,logk[2]), pow(10,logk[3])      
        # tansition rates: G1->S, S->G2, G2->M, (G1, S, G2, M)->D
        dxdt =  [-(kG1S + kphi)*G1 + 2.0*kG2MG1*G2_plus_M,
                -(kSG2M + kphi)*S + kG1S*G1,
                -(kG2MG1 + kphi)*G2_plus_M + kSG2M*S,
                (G1 + S + G2_plus_M)*kphi]
        return dxdt
    return  odeint(cell_cycle_noM, x0, t, args=(logk,), 
                   rtol=1e-4, atol=1e-4, hmin=1e-4)


def obj_func(logrates, dfdata, phases=['G1', 'S', 'G2', 'M', 'D'], func='sum_sq'):
    '''
    Objective function to be minimized.
    Pararmeters:
    -----------

    logrates : list, dtype=float
        rate parameters in log-space.
    dfdata : pandas dataframe
        input data
    func : str
        type of objective function, default = 'sum_sq'
    '''
    
    x0 = dfdata[dfdata.timepoint==0][phases].values.tolist()[0]
    timepoint = dfdata['timepoint'].tolist()
    if len(phases)==5: 
        xmodel = ode_soln(x0, timepoint, logrates)
    elif len(phases)==4:
        xmodel = ode_soln_noM(x0, timepoint, logrates)
    dfmodel = pd.DataFrame()
    for t in range(len(timepoint)):
        dfmodel = dfmodel.append(pd.DataFrame(np.array([xmodel[t]]), columns=phases))
    dfmodel = dfmodel.reset_index().drop(columns='index')
    dfmodel['timepoint'] = timepoint
    if func=='sum_sq':
        sumsq = 0
        s = (dfmodel-dfdata)[phases].values.flatten()
        for s in s: 
            sumsq+=s*s
        return sumsq
    

def fit_rates(dfdata,
              phases=['G1', 'S', 'G2', 'M', 'D'],
              lower_bound=np.log10(10e-8),
              upper_bound=np.log10(10e1),
              num_itr=100, 
              max_obj=np.inf,
              method=None):
    '''
    Fits the phase transistion rates for a given set of time-points and corresponding counts in cell cycle phases 
    
    Parameters:
    -----------
    dfdata: pandas.Dataframe
        Dataframe consisting of counts in cell cycle phases
    phases: array, str
        List of the cell cycle phases included in the model.
        default=['G1', 'S', 'G2_plus_M', 'D'] 
    lower_bound: float
        Lower bound (in log-space) on the rates, default=10e-8
    upper_bound: float
        Upper bound (in logspace) on the rates, default=10e1        
    num_itr: int
        Number of iterations, independent initial guesses for scipy.minimize(), default=100
    max_obj: float
        Maximum objective function
    method: str
        Specifies the method used in spcipy.minimize()
    
    Returns:
    --------
    k: array, float
        the rate parameters 
    '''
    bounds = [tuple([lower_bound, upper_bound]) for i in range(len(phases))]
    bounds = tuple(bounds)

    for itr in range(num_itr):
        #print('Iteration = %s' % itr)
        seed(uniform(0,100000))
        logk0 = [uniform(lower_bound, upper_bound) for ph in range(len(phases))]
        try:    
            rates = minimize(obj_func,
                            logk0,
                            args=(dfdata, phases), 
                            method=method, 
                            bounds=bounds)
        except OverflowError:
            print ('OverflowError')
            print ('Iteration = %s' % itr)
        if rates.fun < max_obj:
            max_obj = rates.fun
            k = [pow(10, y) for y in rates.x.tolist()]
    return (k, max_obj)


def get_rates(cell,
              drug,
              input_file='INPUT/cell_cycle_data.csv',
              num_itr=100,
              phases = ['G1', 'S', 'G2', 'M', 'D'],
              params=[r'$k_{S}$', r'$k_{G2}$',r'$k_{M}$', r'$k_{G1}$', r'$k_{\phi}$']
            ):
    '''
    Estimates the cell cycle phase transistion rates based on the deep dye drop for a given cell-drug combination.

    Parameters:
    -----------
    cell: str
        cell-line
    drug: str
        drug
    input_file: str
        Input file consisting of a long table of cell, drug, concentrations, countsin cell cycle phases.
    num_itr: int
        Number of iterations, idependent initial guess for scipy.minimize()
    phases: list, str
        The cell cycle phases for which the rates need to be estimated
    params: list, str
        The labels for the rate parameters
    
    Returns:
    --------
    dfout: pandas.DataFrame
        DataFrame consisting of the fitted rates at various concentrations for the given cell-drug combination 
    '''
    df = pd.read_csv(input_file)
    df['concentration'] = [round(c, 6) for c in df.concentration.tolist() ]
    assert not df[(df.agent==drug) & (df.cell_line==cell)].empty, \
                'Cell-drug combination not present in the dataset.'
    input_columns = ['timepoint', 'concentration', 'well'] + phases
    output_columns = ['cell_line', 'agent', 'concentration', 
                      'parameter', 'rate', 'obj_func']
    dfctrl = df[(df.cell_line==cell) & 
                (df.agent=='DMSO') & 
                (df.timepoint==0)][input_columns]
    dftreat = df[(df.cell_line==cell) & (df.agent==drug)][input_columns]
    concentration = dftreat.concentration.unique().tolist()
    grpconc = dftreat.groupby('concentration', as_index=False)    
    dfout = pd.DataFrame(columns=output_columns)
    for conc in concentration:
        print ('%s %s = %2.4f uM'%(cell, drug.split('_')[0], conc))
        dfdata = pd.DataFrame()
        dfdata = dfdata.append(dfctrl)
        dfdata = dfdata.append(grpconc.get_group(conc))
        dfdata = dfdata.reset_index().drop(columns=['index', 'concentration', 'well'])
        dfdata = dfdata.groupby('timepoint', as_index=False).mean()
        k, max_obj = fit_rates(dfdata, phases=phases, num_itr=num_itr)
        for i in range(len(params)):
            dfout = dfout.append(pd.DataFrame([[cell,
                                                drug,
                                                conc,
                                                params[i],
                                                k[i],
                                                max_obj]],
                                                columns=output_columns))
    return dfout


def get_rates_rep(cell,
              drug,
              input_file='INPUT/cell_cycle_data.csv',
              num_itr=100,
              phases = ['G1', 'S', 'G2', 'M', 'D'],
              params=[r'$k_{S}$', r'$k_{G2}$', r'$k_{M}$' r'$k_{G1}$', r'$k_{\phi}$']
            ):
    df = pd.read_csv(input_file)
    df['concentration'] = [round(c, 6) for c in df.concentration.tolist() ]
    assert not df[(df.agent==drug) & (df.cell_line==cell)].empty, \
                'Cell-drug combination not present in the dataset.'
    input_columns = ['timepoint', 'concentration', 'well'] + phases
    output_columns = ['cell_line', 'agent', 'concentration',
                      'parameter', 'rate', 'obj_func']
    dfctrl = df[(df.cell_line==cell) &
                (df.agent=='DMSO') &
                (df.timepoint==0)][input_columns]
    dftreat = df[(df.cell_line==cell) & (df.agent==drug)][input_columns]
    concentration = dftreat.concentration.unique().tolist()
    grpconc = dftreat.groupby('concentration', as_index=False)
    dfout = pd.DataFrame(columns=output_columns)
    for conc in concentration:
        print ('%s %s = %2.4f uM'%(cell, drug.split('_')[0], conc))
        dfconc = grpconc.get_group(conc)
        for r in range(dfconc.shape[0]):
            dfdata = pd.DataFrame()
            dfdata = dfdata.append(dfctrl)
            dfdata = dfdata.append(dfconc.iloc[r])
            dfdata = dfdata.reset_index().drop(columns=['index',
                                                        'concentration', 
                                                        'well']
                                               )
            dfdata = dfdata.groupby('timepoint', as_index=False).mean()
            print (dfdata)
            if dfdata.dropna().shape[0]>1:
                k, max_obj = fit_rates(dfdata, phases=phases, num_itr=num_itr)
                for i in range(len(params)):
                    dfout = dfout.append(pd.DataFrame([[cell,
                                                        drug,
                                                        conc,
                                                        params[i],
                                                        k[i],
                                                        max_obj]],
                                                        columns=output_columns)
                                        )
            else: continue
    return dfout


def plot_rates(df, out_dir=None):
    sns.set_style('darkgrid')
    sns.set_context('talk', font_scale=0.75)
    ax = sns.lineplot(data=df,
                      x='concentration', 
                      y='rate', 
                      hue='parameter', 
                    )
    ax.set_xscale('log')
    ax.set_xlabel(r'log$_{10}$C [$\mu$M]')
    ax.set_ylabel(r'Rate [hr$^{-1}$]')
    plt.title('%s, %s' %(df.cell_line.unique()[0], df.agent.unique()[0].split('_')[0]))
    plt.legend(bbox_to_anchor=(1.0, 0.75), loc=2, fontsize=10, frameon=False)
    plt.tight_layout()
    if out_dir is not None:
        plt.savefig('%s%s_%s.pdf' % (out_dir, 
                                    df.cell_line.unique()[0], 
                                    df.agent.unique()[0]
                                    )
                    )
        plt.close()
    return None
