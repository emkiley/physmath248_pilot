'''
This module provides a simple Gaussian elimination routine. 

This solves simple sets of linear equations. More advanced 
aspects, like pivoting, were not covered in class and are 
therefore not implemented.
'''

from numpy import *

def dot_prod(u,v):
    '''
    vector product
    
    input: two vectors of equal lenght
    output: scalar 
    '''
    if not len(u) == len(v): 
        print('Error in matrix.dot_prod: vectors have not the same length')
        return 
    uv_prod = 0
    for i in range(len(u)):
        uv_prod += u[i]*v[i]
    return uv_prod

def mat_vec_prod(A,v):
    '''
    returns matrix-vector product A*v
    
    Parameters:
    -----------
    A : array of arrays
       matrix
    v : array
    
    Returns:
    --------
        array, float
        vector A*v
        
    Notes:
    ------
    https://en.wikipedia.org/wiki/Matrix_multiplication    
    '''
    
    # tests
    # len(v) == number of columns of A, which is len(A.T)    
    if not len(v) == len(A.T): 
        print('Error in matrix.mat_vec_prod: not len(vector) == len(A.T)')
        return 

    # output vector has length = numper of rows of A
    u = zeros(len(A))
    for i in range(len(u)):
        # the ith component of u is the dot product of the ith row of A
        # with the vector v
        u[i] = dot_prod(A[i],v)
    return u

def gauss_elim(A,u):
    '''
    solves equation A*v = u by Gaussian elimination
    
    Parameters:
    -----------
    A : array of arrays
       coefficient matrix
    u : array
        RHS vector of set of equations
        len(u) == number of columns in A
    Returns:
    --------
    v : array, float
        vector v
    BB : array of arrays
       coefficient matrix with RHS in diagonal form
        
    Notes:
    ------
    A must be square matrix.    
    '''
    
    # tests
    # len(u) == number of rows of A, which is len(A.T)    
    if not len(u) == len(A.T): 
        print('Error in matrix.gaus_elim: not len(vector) == len(A.T)')
        return 
    # setup work copy of A, dtype=float and add RHS
    AA = vstack((A.T,u)).T
    if not 'float' in str(AA.dtype):
        AA=AA.astype(float) 
    # iter over rows from first (index 0) to second to last (index len(A)-1)
    for i in range(len(A)):
        AA[i] = AA[i]/AA[i,i]
        for j in range(i+1,len(AA)):
            AA[j] -= AA[j,i]*AA[i]
    # this should give us the diagonal form 
    # return AA
    BB = copy(AA)
    
    # back-substitution
    m  = len(A)-1     # highest row/col index (we do only square matrices)
    v  = zeros(m+1,float) 
    u  = AA.T[-1]       # extract RHS
    AA = delete(AA,m+1,1) # recover diagonalized coefficient matrix  
    for j in range(m,-1,-1):
        v[j] = u[j] - dot_prod(v,AA[j])
    return BB,v
