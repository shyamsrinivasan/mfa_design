# script to test Costas's MFA design problem
# problem is MILP where fluxes are known, choices of mutants are integer decision variables
from gurobipy import *
import numpy as np
import pandas as pd
import os.path
from input.custom_functions import calc_sse, read_data

# read data from excel file
# file_name = os.path.join(os.getcwd(), 'e_coli_core_flux.xlsx')
# read data from csv file
file_name = os.path.join(os.getcwd(), './output/pfba_ecoli_core_filtered.csv')
fluxes, _ = read_data(file_name, file_type='csv')
# rxn_id = flux_df['Abbreviation'].values.tolist()

# fluxes = [np.array([1, .5, .5, .5, .5, .5, .5]), np.array([1, 1, 0, 1, 0, 1, 0]), np.array([1, 0, 1, 0, 1, 0, 1])]
n_mutants = len(fluxes)
# number of mutants (# experiments required)
req_n_mutants = 2
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
# m.write("./output/example.lp")

# set initial value of first two mutants to 1
mut_select[0].start = 1
mut_select[1].start = 1
mut_mult[0].start = 1.0

# optimize model
m.optimize()

nnz_mutant = []
print('\nFinal Solution:\n')
for i_mutant in range(n_mutants):
    # print all solutions to screen
    print('Mutant {} is included {}\n'.format(i_mutant, mut_select[i_mutant].x))
    # collect all non-zero mutant indices
    if mut_select[i_mutant].x > 0:
        nnz_mutant.append(i_mutant)

# collect combos containing nnz_mutants
considered_combos = [list(map(lambda x: True if i_mutant == x[1] or i_mutant == x[2] else False, cont_index)) for i_mutant in nnz_mutant]
# actual_combos = []
# for i_mutant, i_mut_info in enumerate(considered_combos):
#     for j_combo in cont_index:
#         if list(i_mut_info)[j_combo[0]] and mut_mult[j_combo[0]].x:
#             print('Mutant ID: {} Combo ID: {}\n'.format(nnz_mutant[i_mutant], j_combo[0]))
#             actual_combos.append({'combo_id': j_combo[0], 'mutant_ids': j_combo[1:],
#                                   'mip_value': mut_mult[j_combo[0]].x})
actual_combos = [{'combo_id': j_combo[0], 'mutant_ids': j_combo[1:], 'mip_value': mut_mult[j_combo[0]].x}
                 for i_mutant, i_mut_info in enumerate(considered_combos) for j_combo in cont_index
                 if list(i_mut_info)[j_combo[0]] and mut_mult[j_combo[0]].x]
# print all continuous variables - debug only
print('\n Combinations in Solution:\n')
for i_combo in cont_index:
    print('Combination {} - Mutant 1 = {} and Mutant 2 = {} - Included: {}\n'.format(i_combo[0], i_combo[1],
                                                                                     i_combo[2],
                                                                                     mut_mult[i_combo[0]].x))

print('\nNone\n')
