import numpy as np
import matplotlib.pyplot as plt
import math

#%% Tridiagonal Matrix Solver

def tridiag_solver(A, d):
    
    # create the arrays for each diagonal
    a = np.array(np.diagonal(A, offset = -1))
    b = np.array(np.diagonal(A))
    c = np.array(np.diagonal(A, offset = 1))
    
    # find the length of the solution vector
    n = np.prod(d.shape)
    
    # modify the coefficients
    for i in range(1, n):
        m = a[i-1]/b[i-1]
        b[i] = b[i] - m*c[i-1]
        d[i] = d[i] - m*d[i-1]
    
    # create the solution vector and assign x[-1]
    x = np.array(b)
    x[-1] = d[-1]/b[-1]
    
    # loop through the x vector to find the approximations
    for k in range(n-2, -1, -1): # start at n-2, increment by -1
        x[k] = (d[k]-c[k]*x[k+1])/b[k]
    
    # return the solution vector
    return x

# test tridiag_solver function
A = np.array([[10,2,0,0],[3,10,4,0],[0,1,7,5],[0,0,3,4]],dtype=float)
dd = np.array([3,4,5,6.])
print('Tridiagonal matrix:')
print('A =', A)
print('Right-hand side vector:')
print('d =', dd)

x = tridiag_solver(A, dd)
print('Approximation to the solution vector:')
print('x =', x)
print('Product of matrix A and approximation x:')
print('d =', A.dot(x))

#%% Differential Equations Solver - developed from tridiag_solver function above

def diff_eq(f, aa, bb, alpha, beta, N):
    # input variables: f (differential equation), aa and bb (endpoints), alpha and beta (boundary conditions), N (number of unknowns)
    
    h = (bb-aa)/N # calculate the step size
    xx = np.linspace(aa,bb,N) # pre-allocate the x values
    d = np.array(f(xx)) # calculate the right-hand side vector
    
    # create the tridiagonal matrix
    A = np.zeros([N,N])
    for i in range(N):
        for k in range(N):
            if i==k:
                A[i,k] = -2/h**2
            elif i==k-1:
                A[i,k] = 1/h**2
            elif i==k+1:
                A[i,k] = 1/h**2
    
    # create the arrays for each diagonal of the matrix
    a = np.array(np.diagonal(A, offset = -1))
    b = np.array(np.diagonal(A))
    c = np.array(np.diagonal(A, offset = 1))
    
    # find the length of the solution vector
    n = np.prod(d.shape)
    
    # modify the coefficients
    for i in range(1, n):
        m = a[i-1]/b[i-1]
        b[i] = b[i] - m*c[i-1]
        d[i] = d[i] - m*d[i-1]
    
    # create the solution vector and use the boundary condition to find the last element of the vector
    u = np.zeros(N)
    u[-1] = beta
    
    for k in range(n-2, -1, -1):
        if k==0:
            # use the boundary condition for the first element of the vector
            u[k] = alpha
        else:
            u[k] = (d[k]-c[k]*u[k+1])/b[k]
    
    # return the solution vector
    return u

# Trial run:
# f = lambda x: np.sin(x)
# u = lambda x: -np.sin(x)
# uu = diff_eq(f, u, 0, math.pi, 15)

# Example one:
# f = lambda x: 20*(x**3)
# u = lambda x: x**5
# N = 30
# aa = 0
# bb = 2
# alpha = u(aa)
# beta = u(bb)

# Example two:  
f = lambda x: -np.sin(x) + 6*x
u = lambda x: np.sin(x) + x**3
N = 50
aa = (math.pi)/4
bb = 2*math.pi
alpha = u(aa)
beta = u(bb)

uu = diff_eq(f, aa, bb, alpha, beta, N)

xx = np.linspace(aa, bb, N)

error = abs(u(xx)-uu)
print('Approximation to the solution:')
print('u =', uu)
print('Error values at each approximation:')
print('Error =', error)

plt.plot(xx,uu,label='numerical solution')
plt.plot(xx,u(xx), label='exact solution')

plt.xlabel('x')
plt.ylabel('y')

plt.title('Solution by Tridiagonal Algorithm')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.show()

#%% Pentadiagonal Matrix Solver

def pentadiag_solver(A, v):
    # create the arrays for each diagonal
    g = np.concatenate((np.array([0,0]),np.array(np.diagonal(A, offset = -2))))
    h = np.concatenate((np.array([0]),np.array(np.diagonal(A, offset = -1))))
    d = np.array(np.diagonal(A, offset = 0))
    e = np.concatenate((np.array(np.diagonal(A, offset = 1)),np.array([0])))
    f = np.concatenate((np.array(np.diagonal(A, offset = 2)),np.array([0,0])))
    
    m = np.size(A,1)
    n = m-1
    
    # pre-allocate the vectors that will be used
    alpha = np.zeros(m)
    gamma = np.zeros(m-1)
    delta = np.zeros(m-2)
    beta = np.zeros(m)
    
    c = np.zeros(m)
    x = np.zeros(m)
    
    # assign the elements to the vectors to get A = LR
    alpha[0] = d[0]
    gamma[0] = e[0]/alpha[0]
    delta[0] = f[0]/alpha[0]
    beta[1] = h[1]
    alpha[1] = d[1] - beta[1]*gamma[0]
    gamma[1] = (e[1]-beta[1]*delta[0])/alpha[1]
    delta[1] = f[1]/alpha[1]
    
    for k in range(2,n-1):
        beta[k] = h[k]-g[k]*gamma[k-2]
        alpha[k] = d[k]-g[k]*delta[k-2]-beta[k]*gamma[k-1]
        gamma[k] = (e[k]-beta[k]*delta[k-1])/alpha[k]
        delta[k] = f[k]/alpha[k]
        
    beta[n-1] = h[n-1]-g[n-1]*gamma[n-3]
    alpha[n-1] = d[n-1] - g[n-1]*delta[n-3]-beta[n-1]*gamma[n-2]
    gamma[n-1] = (e[n-1]-beta[n-1]*delta[n-2])/alpha[n-1]
    beta[n] = h[n]-g[n]*gamma[n-2]
    alpha[n] = d[n]-g[n]*delta[n-2]-beta[n]*gamma[n-1]
    
    # find the c vector, where v = Lc
    c[0] = v[0]/alpha[0]
    c[1] = (v[1]-beta[1]*c[0])/alpha[1]
    
    for j in range(2,m):
        c[j] = (v[j]-g[j]*c[j-2]-beta[j]*c[j-1])/alpha[j]
    
    # use backwards substitution to find the x vector, where c = Rx
    x[n] = c[n]
    x[n-1] = c[n-1]-gamma[n-1]*x[n]
    
    for i in range(n-2, -1, -1):
        x[i] = c[i]-gamma[i]*x[i+1]-delta[i]*x[i+2]
    
    return x

A = np.array([[7,3,2,0,0,0],[3,6,7,5,0,0],[2,7,9,4,6,0],[0,5,4,8,9,5],[0,0,6,9,3,4],[0,0,0,5,4,2]])
v = np.array([23,70,95,99,83,36])

print('Pentadiagonal matrix: A =')
print(A)
print('Right-hand side vector:')
print('d =', v)

x = pentadiag_solver(A, v)
print('Approximation to the solution vector:')
print('x =', x)
print('Product of matrix A and approximation x:')
print('d =', A.dot(x))

error = abs(v-A.dot(x))
print('Error values at each approximation:')
print('error:', error)
print('Maximum error:', error.max())
