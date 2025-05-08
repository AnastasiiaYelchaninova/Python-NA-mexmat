import numpy as np

# solves a system of linear algebraic equations (SLAE) using Gauss-Seidel algorithm

A = np.array([[5.,3.,0.,1.],
           [2.,-4.,1.,0.],
           [4.,3.,12.,-1.],
           [3.,0.,1.,6.]])
b = np.array([[1.],
           [2.],
           [6.],
           [4.]])

print("Coefficients matrix:")
print(A)
print("Right side vector:")
print(b)
print()

tolerance = 1e-6
max_iterations = 50

def gauss_seidel(A, b, tolerance, max_iterations):
    x = np.zeros_like(b, dtype=np.double)
    for k in range(max_iterations):
        x_old  = x.copy()
        for i in range(A.shape[0]): # Gauss-Seidel algorithm: calculate a new approximate value x[i]
            x[i] = (b[i] - np.dot(A[i,:i], x[:i]) - np.dot(A[i,(i+1):], x_old[(i+1):])) / A[i ,i]
        if np.linalg.norm(x - x_old, ord=np.inf) / np.linalg.norm(x, ord=np.inf) < tolerance: # break conditions 
            print("Iterations:", k)
            print()
            break        
    return x

x = gauss_seidel(A, b, tolerance, max_iterations)
print(x)
print()

nev = np.dot(A,x)-b
print("Residual vector:")
print(nev)
