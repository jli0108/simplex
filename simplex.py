import numpy as np
import sys
# -- Solve your LP problems with this simple code!!! -----
# -- Your LP must be in the following format: ------------
# -- Maximize/Minimize c^T x -----------------------------
# -- subject to Ax <= b ----------------------------------

# -- Modify if you want to maximize or minimize ----------
maximize = True
# -- Modify these arrays in the correct format -----------
c = np.array([[10, -57, -9, -24]])

A_B = np.array([[0.5, -5.5, -2.5, 9],
                [0.5, -1.5, -0.5, 1],
                [1, 0, 0, 0]])

b = np.array([[0],
              [0],
              [1]])

assert A_B.shape[0] == b.shape[0]
assert c.shape[1] == A_B.shape[1]
# solves problem in canonical form
def solve_canonical(c, A_B, b, maximize):
    if not maximize:
        c = -c
    tableau, basis_variables = maximize_canonical(c, A_B, b)
    print("----- Solution -----")
    if maximize:
        print("Max value of objective function:", -tableau[0,-1])
    else:
        print("Min value of objective function:", tableau[0,-1])
    for i in range(basis_variables.shape[0]):
        print(f"x{basis_variables[i]+1} = {tableau[i+1,-1]}")
    print("Set all other variables to zero.")

# solves maximization problem in canonical form
# returns the final tableau and a list of the basis variables
def maximize_canonical(c, A_B, b):
    if b.min() < 0:
        sys.exit("Initial feasibility problem. Not implemented yet.")

    m, n = A_B.shape
    A_N = np.eye(m)

    A = np.concatenate((A_B, A_N), axis=1)
    c = np.concatenate((c, np.zeros((1,m))), axis=1)
    basis_variables = np.arange(n, m+n)
    tableau = np.concatenate((np.concatenate((c, np.array([[0]])), axis=1), np.concatenate((A, b), axis=1)))

    # choose entering variable
    entering_variable = None
    i = 0
    while i < m+n :#and entering_variable is None:
        # Uncomment lines above and below to use Bland's rule 
        # (if so, we do not need to check if entering_variable is None)
        #if tableau[0,i] > 0:
        # It should usually be faster if we take the largest coefficient
        if tableau[0,i] > 0 and (entering_variable is None or tableau[0,i] > tableau[0,entering_variable]):
            entering_variable = i
        i = i + 1
    print(tableau)
    while entering_variable is not None:
        
        index_of_leaving_variable = None
        for i in range(tableau[1:m+1,-1].size):
            if (tableau[i+1,entering_variable] > 0 and # only consider if coefficient is positive
                (tableau[i+1,-1] / tableau[i+1,entering_variable] >= 0) and # bound is nonnegative
                (index_of_leaving_variable is None or 
                # find smallest nonnegative coefficient to determine leaving variable (using multiplication for performance)
                tableau[index_of_leaving_variable+1,entering_variable] * tableau[i+1,-1] < tableau[i+1,entering_variable] * tableau[index_of_leaving_variable+1,-1] )):
                index_of_leaving_variable = i
        
        print(f"x{entering_variable+1} enters, x{basis_variables[index_of_leaving_variable]+1} leaves")
        # variable enters the basis
        basis_variables[index_of_leaving_variable] = entering_variable

        # pivot operation
        tableau[index_of_leaving_variable+1,:] = tableau[index_of_leaving_variable+1,:] / tableau[index_of_leaving_variable+1,entering_variable]

        for i in range(tableau.shape[0]):
            if i != index_of_leaving_variable + 1:
                tableau[i,:] -= tableau[i,entering_variable] * tableau[index_of_leaving_variable+1,:]
        
        print(tableau)
        # choose entering variable
        entering_variable = None
        i = 0
        # use Bland's rule
        while i < m+n and entering_variable is None:
            if tableau[0,i] > 0:
            #if tableau[0,i] > 0 and (entering_variable is None or tableau[0,i] > tableau[0,entering_variable]):
                entering_variable = i
            i = i + 1
    return tableau, basis_variables

solve_canonical(c, A_B, b, maximize)