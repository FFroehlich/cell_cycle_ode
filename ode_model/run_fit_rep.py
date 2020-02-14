import sys
import pandas as pd
from ode_model import get_rates_rep, plot_rates

drug = sys.argv[1]
input_file = sys.argv[2]
output_dir = sys.argv[3]

df = pd.read_csv(input_file)
dfout = pd.DataFrame()

for cell in df[df.agent==drug]['cell_line'].unique().tolist():
    print (drug, cell)
    dffit = get_rates_rep(cell, drug, input_file, num_itr=100)
    dffit.to_csv('%s%s_%s.csv'%(output_dir, cell, drug), index=False)
    #plot_rates(dffit[dffit.obj_func<10e3], out_dir=output_dir)
    dfout = pd.concat([dfout, dffit], sort=False, ignore_index=True)
    #dfout = dfout.append(dffit)
dfout = dfout.reset_index()
dfout = dfout.drop(columns='index')
dfout.to_csv('%s%s.csv' % (output_dir, drug), index=False)

