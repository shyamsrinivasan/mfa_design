# script to test Costas's MFA design problem
# problem is MILP where fluxes are known, choices of mutants are integer decision variables
from gurobipy import *
import numpy as np
import pandas as pd
import os.path


def calc_sse(val_1, val_2):
    # calculate lsq distance between val_1 and val_2
    sse_estimate = np.sqrt(np.sum((abs(val_1) - abs(val_2))**2)/len(val_1))
    return sse_estimate


def read_data(file_name, file_type='excel'):

    df = []
    if file_type == 'excel':
        # read data from excel file
        df = pd.read_excel(file_name, sheet_name='Sheet4')
    elif file_type == 'csv':
        # read data from csv file
        df = pd.read_csv(file_name, sep='\t')  # , lineterminator='\r\n')

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


# read data from excel file
# file_name = os.path.join(os.getcwd(), 'e_coli_core_flux.xlsx')
# read data from csv file
file_name = os.path.join(os.getcwd(), 'pfba_ecoli_core.csv')
fluxes, _ = read_data(file_name, file_type='csv')
# rxn_id = flux_df['Abbreviation'].values.tolist()

# fluxes = [np.array([1, .5, .5, .5, .5, .5, .5]), np.array([1, 1, 0, 1, 0, 1, 0]), np.array([1, 0, 1, 0, 1, 0, 1])]
n_mutants = len(fluxes)
# number of mutants (# experiments required)
req_n_mutants = 5
# pre-calculate sse estimates for different flux pairs
flux_sse = [calc_sse(fluxes[i], fluxes[j]) for i in range(0, n_mutants-1) for j in range(i, n_mutants) if i != j]
combo_index = [(i, j) for i in range(0, n_mutants-1) for j in range(i, n_mutants) if i != j]
cont_index = [(combo_ind, i_combo_index[0], i_combo_index[1])for combo_ind, i_combo_index in enumerate(combo_index)]

# create empty model
m = Model("mfa_design")

# variables
# binary variable(s) - choice of mutant, bounds: y is binary (0, 1)
mut_select = m.addVars(range(n_mutants), vtype=GRB.BINARY, name='mut_select')
# continuous variable(s) - product of binary variables for use in objective function, bounds: 0 <= r <= 1
mut_mult = m.addVars(range(len(flux_sse)), lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name='mut_select_cont')

# objective: (vi - vj)**2 * rij i-> 1 to K-1 and j-> 2 to K
m.setObjective(quicksum(flux_sse[i_flux_comb] * mut_mult[i_flux_comb] for i_flux_comb in range(len(flux_sse))),
               sense=GRB.MAXIMIZE)

# constraint 1: rij <= yi
test_c1 = [(j_cont_index[0], j_cont_index[1]) for j_cont_index in cont_index]

m.addConstrs((mut_mult[j_cont_index[0]] - mut_select[j_cont_index[1]] <= 0 for j_cont_index in cont_index), "c1")
# constraint 2: rij <= yj
m.addConstrs((mut_mult[j_cont_index[0]] - mut_select[j_cont_index[2]] <= 0 for j_cont_index in cont_index), "c2")
# constraint 3: rij >= yi + yj - 1
m.addConstrs((mut_select[j_cont_index[1]] + mut_select[j_cont_index[2]] - mut_mult[j_cont_index[0]] - 1 <= 0
              for j_cont_index in cont_index), "c3")
# constraint 4: sum(y) = L (L is fixed)
m.addConstr(quicksum(mut_select) - req_n_mutants, rhs=0, name="c4", sense=GRB.LESS_EQUAL)

# force model update
m.update()

# Save model to file
m.write("example.lp")

# optimize model
m.optimize()

print('\nNone\n')