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
c : np.ndarray = np.array([[7, 6, 5, -2, 3]])

A_B : np.ndarray = np.array([[1, 3, 5, -2, 2],
                [4, 2, -2, 1, 1],
                [2, 4, 4, -2, 5],
                [3, 1, 2, -1, -2]])

b : np.ndarray = np.array([[4],
              [3],
              [5],
              [1]])

assert A_B.shape[0] == b.shape[0]
assert c.shape[1] == A_B.shape[1]

def solve_LP(c : np.ndarray, A_B : np.ndarray, b : np.ndarray, maximize : bool, standard : bool) -> None:
    if (standard):
        solve_standard(c, A_B, b, maximize)
    else:
        solve_canonical(c, A_B, b, maximize)

# solves problem in canonical form
def solve_canonical(c : np.ndarray, A_B : np.ndarray, b : np.ndarray, maximize : bool) -> None:
    if not maximize:
        c = -c
    optimal_cost, RHS, basic_variables, shadow_prices, A_B_inv = maximize_canonical(c, A_B, b)
    print("----- Solution -----")
    if maximize:
        print("Max value of objective function:", optimal_cost)
    else:
        print("Min value of objective function:", -optimal_cost)
    for i in range(basic_variables.shape[0]):
        print(f"x{basic_variables[i]+1} = {RHS[i,0] : 0.7f}")
    for i in range(shadow_prices.shape[1]):
        print(f"Shadow price of x{basic_variables[i]+1}: {shadow_prices[0,i] : 0.7f}")
    for i in range(RHS.shape[0]):
        allowable_increase = None
        allowable_decrease = None
        for j in range(RHS.shape[0]):
            if A_B_inv[j,i] > 0:
                if allowable_decrease is None or allowable_decrease > RHS[i,0] / A_B_inv[j,i]:
                    allowable_decrease = RHS[i,0] / A_B_inv[j,i]
            elif A_B_inv[j,i] < 0:
                if allowable_increase is None or allowable_increase > - RHS[i,0] / A_B_inv[j,i]:
                    allowable_increase = - RHS[i,0] / A_B_inv[j,i]
        if allowable_decrease is None:
            allowable_decrease = np.Infinity
        if allowable_increase is None:
            allowable_increase = np.Infinity      
        print(f"Row {i+1}: Allowable decrease = {allowable_decrease : 0.001f}, Allowable increase = {allowable_increase : 0.001f}")
    #print(A_B_inv)
    print("Set all other variables to zero.")

def solve_standard(c, A_B, b, maximize):
    sys.exit("Not implemented without tableaus yet. Use simplex.py")

# takes problem in standard form and solves phase 1 LP
def solve_phase_one(A_B, b):
    sys.exit("Not implemented without tableaus yet. Use simplex.py")

# solves maximization problem in canonical form
# returns the final tableau and a list of the basis variables
def maximize_canonical(c_N, A_N, b):
    #if b.min() < 0:
    #    sys.exit("Initial feasibility problem. Not implemented yet.")

    m, n = A_N.shape
    A_B = np.eye(m)
    for i in range(m):
        if b[i] < 0:
            A_N[i] *= -1
            b[i,0] *= -1
            A_B[i] *= -1
    A = np.concatenate((A_N, np.eye(m)), axis=1)
    #print(A)
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
    while i < n:# and index_of_entering_variable is None:
        # Uncomment lines above and below to use Bland's rule 
        # (if so, we do not need to check if entering_variable is None)
        if reduced_c_N[0,i] > 0 and (index_of_entering_variable is None or reduced_c_N[0,i] > reduced_c_N[0,index_of_entering_variable]):
        # It should usually be faster if we take the largest coefficient
        #if c_N[i] > 0 and (entering_variable is None or c_N[i] > c_N[entering_variable]):
            index_of_entering_variable = i
        i = i + 1

    while index_of_entering_variable is not None:

        reduced_col_of_A_N = np.dot(A_B_inv, A_N[:,index_of_entering_variable])

        #print(A_B_inv)
        #print(b_bar[:,0] / reduced_col_of_A_N)
        index_of_leaving_variable = None
        for i in range(m):
            if (reduced_col_of_A_N[i] > 0 and
                (index_of_leaving_variable is None or
                b_bar[i,0] * reduced_col_of_A_N[index_of_leaving_variable] < b_bar[index_of_leaving_variable,0] * reduced_col_of_A_N[i])):

                index_of_leaving_variable = i

        if index_of_leaving_variable is None:
            sys.exit("LP is unbounded")
        entering_variable = nonbasic_variables[index_of_entering_variable]
        leaving_variable = basic_variables[index_of_leaving_variable]
        #print("basis",basic_variables+1)
        #print("nonbasic",nonbasic_variables+1)
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
        #print(A_B)
        #if (np.linalg.det(A_B) == 0):

        A_B_inv = inv(A_B)
        y = np.dot(c_B, A_B_inv)
        b_bar = np.dot(A_B_inv, b)
        reduced_c_N = c_N - np.dot(y, A_N)

        index_of_entering_variable = None
        i = 0
        #print(reduced_c_N)
        # for first iteration, entering variable is equal to index
        while i < n:# and index_of_entering_variable is None:
            # Uncomment lines above and below to use Bland's rule 
            # (if so, we do not need to check if entering_variable is None)
            if reduced_c_N[0,i] > 0 and (index_of_entering_variable is None or reduced_c_N[0,i] > reduced_c_N[0,index_of_entering_variable]):
            # It should usually be faster if we take the largest coefficient
            #if c_N[i] > 0 and (entering_variable is None or c_N[i] > c_N[entering_variable]):
                index_of_entering_variable = i
            i = i + 1
        #print(index_of_leaving_variable)
    return np.dot(y, b)[0,0], b_bar, basic_variables, y, A_B_inv

#solve_canonical(c, A_B, b, maximize)
#solve_standard(c, A_B, b, maximize)
solve_LP(c, A_B, b, maximize, standard)