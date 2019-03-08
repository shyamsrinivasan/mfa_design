# Suggest informative mutants to generate MFA data

Scripts provided in this repository implement a procedure to use stoichiometric metabolic model to find a set(s) of mutants 
that are the most dispersed. In detailed, the sum of pairwise Euclidean distances between mutants' flux vectors are 
maximized for a set of mutants.<br>
<br>
In the first step, mutants (found by knocking out reactions) are filtered, mutants' flux vectors are determined using 
parsimonious FBA (pFBA) (http://msb.embopress.org/content/6/1/390), and Euclidean distances are calculated. The 
following mutants are excluded:<br>
1. Essential reactions determined by essential genes via GPR
2. Exchange and transport reactions (even transport reactions with genes associated)
3. Essential reactions in silico
4. Reactions without associated genes
5. Redundant mutants determined by their flux vectors being identical (sum pairwise < 1e-6)
6. Reactions with zero flux (determined by FVA)

In the second step, MILP formulation is used to determine the optimal set of mutants. Alternative solutions can be found with 
integer cut by setting the "number of solutions" parameter.<br>
<br>
Current implementation with ipython notebook used optlang (https://github.com/biosustain/optlang) which use common syntax 
for common solvers such as GLPK, CPLEX, and Gurobi.
