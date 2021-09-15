import numpy as np

# -- Solve your LP problems with this simple code!!! -----
# -- Your LP must be in the following format: ------------
# -- Maximize c^T x --------------------------------------
# -- subject to Ax <= b ----------------------------------

# -- Modify these arrays in the correct format -----------
c = np.array([[4, 3, 7]])

A_B = np.array([[1, 3, 2],
                [2, 1, 3]])

b = np.array([[120],
              [120]])

m, n = A_B.shape
A_N = np.eye(m)

A = np.concatenate((A_B, A_N), axis=1)
c = np.concatenate((c, np.zeros((1,m))), axis=1)
basis_variables = np.arange(n, m+n)
minus_objective = np.array([[0]])
#print(basis_variables)
#print("A:")
#print(A)
#print("c:", c)
# initial feasible solution of zero.
tableau = np.concatenate((np.concatenate((c, minus_objective), axis=1), np.concatenate((A, b), axis=1)))
#print(tableau)

# choose entering variable
entering_variable = None
for i in range(n+m):
    if tableau[0,i] > 0 and (entering_variable is None or tableau[0,i] > tableau[0,entering_variable]):
        entering_variable = i
print("entering:", entering_variable)
while entering_variable is not None:
    # find smallest nonnegative coefficient to determine leaving variable
    constraints = tableau[1:m+1,-1] / tableau[1:,entering_variable]
    index_of_leaving_variable = None
    for i in range(constraints.size):
        if constraints[i] > 0 and (index_of_leaving_variable is None or constraints[i] < constraints[index_of_leaving_variable]):
            index_of_leaving_variable = i
    #print(constraints)
    #print("leaving:", index_of_leaving_variable)
    basis_variables[index_of_leaving_variable] = entering_variable
    #print(basis_variables)
    tableau[index_of_leaving_variable+1,:] = tableau[index_of_leaving_variable+1,:] / tableau[index_of_leaving_variable+1,entering_variable]

    for i in range(tableau.shape[0]):
        if i != index_of_leaving_variable + 1:
            tableau[i,:] -= tableau[i,entering_variable] * tableau[index_of_leaving_variable+1,:]

    #print(tableau)
    
    entering_variable = None
    for i in range(n+m):
        if tableau[0,i] > 0 and (entering_variable is None or tableau[0,i] > tableau[0,entering_variable]):
            entering_variable = i
    #print("entering:", entering_variable)

#print("final basis variables:", basis_variables)
print("-- Solution --")
print("Max value of objective function:", -tableau[0,-1])
for i in range(basis_variables.shape[0]):
    print(f"x{basis_variables[i]} = {tableau[i+1,-1]}")
print("Set all other variables to zero")
