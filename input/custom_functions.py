def calc_sse(val_1, val_2):
    import numpy as np
    # calculate lsq distance between val_1 and val_2
    sse_estimate = np.sqrt(np.sum((abs(val_1) - abs(val_2))**2)/len(val_1))
    return sse_estimate


def read_data(file_name, file_type='excel'):
    import pandas as pd
    import numpy as np
    df = []
    if file_type == 'excel':
        # read data from excel file
        df = pd.read_excel(file_name, sheet_name='Sheet4')
    elif file_type == 'csv':
        # read data from csv file
        df = pd.read_csv(file_name, sep='\t', index_col=0)  # , lineterminator='\r\n')

    # non-mutant column names (excel): Abbreviation, Description, Reaction, GPR, Lower bound, Upper bound,
    # Objective, Subsystem
    # parse data
    all_column_names = df.columns.values
    non_mutant_column_names = ['Abbreviation', 'Description', 'Reaction', 'GPR', 'Lower bound', 'Upper bound',
                               'Objective', 'Subsystem']
    mutant_column_names = [i_column for i_column in all_column_names if i_column not in non_mutant_column_names]

    # collect mutant flux data
    mutant_flux = [df[i_mutant].values for i_mutant in mutant_column_names]

    # set all flux values below threshold to zero
    adjusted_flux = [np.array(list(map(lambda x: 0 if abs(x) < 1e-6 else x, list(i_mutant_flx))))
                     for i_mutant_flx in mutant_flux]
    return adjusted_flux, df