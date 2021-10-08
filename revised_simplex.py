import numpy as np
from numpy.linalg import inv
import sys
# -- Solve your LP problems with this simple code!!! -----
# -- Your LP must be in the following format: ------------
# -- Maximize/Minimize c^T x -----------------------------
# -- subject to Ax <= b ----------------------------------

# -- Modify if you want to maximize or minimize -----------------------
maximize : bool = True
# -- Modify if your problem is in standard or canonical form ----------
standard : bool = False
# -- Modify these arrays in the correct format ------------------------
c : np.ndarray = np.array([[2, 4, 1, 1]])

A_B : np.ndarray = np.array([[1, 3, 0, 1],
                [2, 1, 0, 0],
                [0, 1, 4, 1]])

b : np.ndarray = np.array([[8],
              [6],
              [6]])

assert A_B.shape[0] == b.shape[0]
assert c.shape[1] == A_B.shape[1]

def solve_LP(c : np.ndarray, A_B : np.ndarray, b : np.ndarray, maximize : bool, standard : bool) -> None:
    if standard:
        solve_standard(c, A_B, b, maximize)
    else:
        solve_canonical(c, A_B, b, maximize)

# solves problem in canonical form
def solve_canonical(c : np.ndarray, A_B : np.ndarray, b : np.ndarray, maximize : bool) -> None:
    if not maximize:
        c = -c
    optimal_cost, RHS, basic_variables, shadow_prices = maximize_canonical(c, A_B, b)
    print("----- Solution -----")
    if maximize:
        print("Max value of objective function:", optimal_cost)
    else:
        print("Min value of objective function:", -optimal_cost)
    for i in range(basic_variables.shape[0]):
        print(f"x{basic_variables[i]+1} = {RHS[i,0]}")
    for i in range(shadow_prices.shape[1]):
        print(f"Shadow price of x{basic_variables[i]+1}: {shadow_prices[0,i]}")
    print("Set all other variables to zero.")

def solve_standard(c, A_B, b, maximize):
    if not maximize:
        c = -c
    tableau, basis_variables = solve_phase_one(A_B, b)
    if tableau[0,-1] == 0:
        m = basis_variables.shape[0]
        n = tableau.shape[1] - m - 1
        
        # should also check all basis variables are not auxiliary...
        for i in range(m):
            assert basis_variables[i] < n, "not implemented yet, but problem is feasible"
        for _ in range(n, n+m):
            tableau = np.delete(tableau, n, 1)
        tableau[0,0:n] = c
        print(tableau)
        print("End of Phase I")
        # should rewrite objective...
        print(basis_variables)
        for i in range(m):
            tableau[0] -= tableau[0,basis_variables[i]] * tableau[i+1]
        print(tableau)
        entering_variable = None
        i = 0
        while i < n :#and entering_variable is None:
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
            while i < n and entering_variable is None:
                if tableau[0,i] > 0:
                #if tableau[0,i] > 0 and (entering_variable is None or tableau[0,i] > tableau[0,entering_variable]):
                    entering_variable = i
                i = i + 1
        print("----- Solution -----")
        if maximize:
            print("Max value of objective function:", -tableau[0,-1])
        else:
            print("Min value of objective function:", tableau[0,-1])
        for i in range(basis_variables.shape[0]):
            print(f"x{basis_variables[i]+1} = {tableau[i+1,-1]}")
        print("Set all other variables to zero.")
    else:
        print("Problem is infeasible.")

# takes problem in standard form and solves phase 1 LP
def solve_phase_one(A_B, b):
    # make RHS nonnegative
    m, n = A_B.shape
    for i in range(m):
        if b[i,0] < 0:
            A_B[i] = -A_B[i]
            b[i]
    
    A_N = np.eye(m)
    A = np.concatenate((A_B, A_N), axis=1)
    c = np.concatenate((np.zeros((1, n)), -np.ones((1, m))), axis=1)
    basis_variables : np.ndarray = np.arange(n, m+n)
    tableau : np.ndarray = np.concatenate((np.concatenate((c, np.array([[0]])), axis=1), np.concatenate((A, b), axis=1)))
    for i in range(1, m+1):
        # write objective in terms of nonbasic variables
        tableau[0] = tableau[0] + tableau[i]

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

# solves maximization problem in canonical form
# returns the final tableau and a list of the basis variables
def maximize_canonical(c_N, A_N, b):
    if b.min() < 0:
        sys.exit("Initial feasibility problem. Not implemented yet.")

    m, n = A_N.shape
    A_B = np.eye(m)
    A = np.concatenate((A_N, np.eye(m)), axis=1)
    # c_B is coefficients from original c (not getting reduced)
    c_B = np.zeros((1,m))
    c = np.concatenate((c_N, c_B), axis=1)
    basic_variables = np.arange(n, m+n)
    nonbasic_variables = np.arange(0, n)

    A_B_inv = np.eye(m)
    y = np.zeros((1,m))
    b_bar = b
    reduced_c_N = c_N
    index_of_entering_variable = None
    i = 0
    # for first iteration, entering variable is equal to index
    while i < n and index_of_entering_variable is None:
        # Uncomment lines above and below to use Bland's rule 
        # (if so, we do not need to check if entering_variable is None)
        if reduced_c_N[0,i] > 0:
        # It should usually be faster if we take the largest coefficient
        #if c_N[i] > 0 and (entering_variable is None or c_N[i] > c_N[entering_variable]):
            index_of_entering_variable = i
        i = i + 1

    while index_of_entering_variable is not None:

        reduced_col_of_A_N = np.dot(A_B_inv, A_N[:,index_of_entering_variable])

        index_of_leaving_variable = None
        for i in range(m):
            if (reduced_col_of_A_N[i] > 0 and
                (index_of_leaving_variable is None or
                b_bar[i,0] * reduced_col_of_A_N[index_of_leaving_variable] < b_bar[index_of_leaving_variable,0]) * reduced_col_of_A_N[i]):

                index_of_leaving_variable = i

        if index_of_leaving_variable is None:
            sys.exit("LP is unbounded")
        entering_variable = nonbasic_variables[index_of_entering_variable]
        leaving_variable = basic_variables[index_of_leaving_variable]
        #print(f"x{entering_variable+1} enters, x{leaving_variable+1} leaves")
        # variable enters the basis
        basic_variables[index_of_leaving_variable] = entering_variable
        nonbasic_variables[index_of_entering_variable] = leaving_variable
        
        # A_B gets updated based on new basic variables
        A_B[:,index_of_leaving_variable] = A[:,entering_variable]
        A_N[:,index_of_entering_variable] = A[:,leaving_variable]
        # c_B gets updated based on new basic variables
        c_B[0,index_of_leaving_variable] = c[0,entering_variable]
        c_N[0,index_of_entering_variable] = c[0,leaving_variable]

        A_B_inv = inv(A_B)
        y = np.dot(c_B, A_B_inv)
        b_bar = np.dot(A_B_inv, b)
        reduced_c_N = c_N - np.dot(y, A_N)

        index_of_entering_variable = None
        i = 0
        # for first iteration, entering variable is equal to index
        while i < n and index_of_entering_variable is None:
            # Uncomment lines above and below to use Bland's rule 
            # (if so, we do not need to check if entering_variable is None)
            if reduced_c_N[0,i] > 0:
            # It should usually be faster if we take the largest coefficient
            #if c_N[i] > 0 and (entering_variable is None or c_N[i] > c_N[entering_variable]):
                index_of_entering_variable = i
            i = i + 1
    
    return np.dot(y, b)[0,0], b_bar, basic_variables, y

#solve_canonical(c, A_B, b, maximize)
#solve_standard(c, A_B, b, maximize)
solve_LP(c, A_B, b, maximize, standard)