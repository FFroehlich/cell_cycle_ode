import  pandas as pd

dfgat = pd.DataFrame()
for chunk in pd.read_csv('../../BrCaLines_profiling/ProfilingPaper/cell_cycle_data/merged/batches1to11_gating.results.csv', chunksize=500):
    df = chunk
    df['D'] = df['cell_count__dead'] + df['corpse_count']       #.div(df['cell_count__total'])
    df[['G1', 'S', 'G2', 'M']] = df[['G1', 'S', 'G2', 'M']].multiply(df['cell_count'], axis=0)
    df['G2_plus_M'] = df['G2'] + df['M']
    dfgat = dfgat.append(df)
dfgat = dfgat.replace('72', 72)
dfgat = dfgat.replace('time0_ctrl', 0)
dfgat = dfgat[['cell_line', 'agent', 'timepoint', 'concentration', 'well', 'G1', 'S', 'G2', 'M', 'G2_plus_M', 'D']]

#Edit namespaces
dfgat['agent'] = [drug.replace('/', '_') for drug in dfgat['agent'].tolist()]
dfgat['agent'] = [drug.replace(' ', '') for drug in dfgat['agent'].tolist()]
dfgat = dfgat.replace({'agent': {'torin2' : 'Torin2'}})
dfgat = dfgat[dfgat.cell_line!='MCF10A']
dfgat['cell_line'] = [cell.upper() for cell in dfgat['cell_line'].tolist()]
dfgat['cell_line'] = [cell.replace('-PR', '_PR') for cell in dfgat['cell_line'].tolist()]
dfgat['cell_line'] = [cell.replace('-HME', '_HME') for cell in dfgat['cell_line'].tolist()]
dfgat['cell_line'] = [cell.replace('-','') for cell in dfgat['cell_line'].tolist()]
dfgat['cell_line'] = [cell.replace(' (GM)', '') for cell in dfgat['cell_line'].tolist()]
dfgat['cell_line'] = [cell.replace('ESVAT', 'EVSAT') for cell in dfgat['cell_line'].tolist()]

dfgat.to_csv('INPUT/cell_cycle_data.csv', index=False)
