from math import*

'''Interpolation by Newton's polynomial of the function
f(x) = x*log(x)+2*x*cos(x)
at equidistant nodes and Chebyshev nodes'''

a = 1
b = 2
n = 20
x = []
f = []
h = (b-a)/n

def function(x):
    return x*log(x)+2*x*cos(x)

def func (x): # calculates function values in nodes
    f = []
    for i in range(n+1):
        f.append(function(x[i]))
    return f

for i in range(n+1):
    x.append(a + i*h)
f = func(x)

'''How shall we calculate the coefficients of Newton's polynomial?
To do this, let's write out a system for the coefficients and solve it using the inverse Gauss method'''

def matrixA(x): # creates a matrix for solving the linear system
    A = []
    d = len(x)
    for j in range(d):
        s = []
        for i in range(d-1):
            s.append(0)
        s.append(1)
        A.append(s)
    for j in range(d-1):
        p = (x[j]-x[d-1])
        for i in range(d-2, j-1, -1):
            A[j][i] = p
            p = p * (x[j]-x[i])
    return A  

def Gauss_back(A, b): # reverse step of the Gaussian method
    d = len(A)
    def scalar_product(A, b, n):
        s = 0
        for i in range(n+1, d):
            s = s + b[i]*A[n][i]
        return s    
    for i in range(d-1, -1, -1):
        b[i] = (b[i] - scalar_product(A, b, i))/A[i][i]
    return b

x.reverse()
f.reverse()
A = matrixA(x) # linear system matrix
g = Gauss_back(A, f[:]) # polynomial coefficients
g.reverse()
x.reverse()
f.reverse()

def polynom_string(g, x): # creating a polynomial
    s = 'f(x) = ' + str(g[0]) + '+'
    p = '*(x-' + str(x[0]) + ')'
    for i in range(1, len(g)):
        s = s + str(g[i]) + p + '+'
        p = p + '*(x-' + str(x[i]) + ')'
    return s[: -1].replace('+-', '-')

def polynom_value(g, x, y): # calculate the value of the polynomial at the given point y
    s = g[0]
    p = (y - x[0])
    for i in range(1, len(g)):
        s = s + g[i]*p
        p = p * (y - x[i])
    return s

pol = [] # array of polynomial values at the partition points
for xi in x:
    pol.append(polynom_value(g, x, xi))

D =[] # array of errors
for i in range(n+1):
    D.append(f[i]-pol[i])

print("f(x) = x[i]*log(x[i])+2*x[i]*cos(x[i])")
print()
print(polynom_string(g ,x))
print()
print("{:^12}{:^12}{:^12}{:^12}".format("x_i", "f_i", "pol_i", "f(x_i)-pol(x_i)"))
for i in range (n+1):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x[i], f[i], pol[i], D[i]))
print()

#all same for "changed" x:
alpha = 0.7
x_1 = []
f_1 = []
pol_1 = []
D_1 =[]

for i in range(n+1):
    x_1.append(x[i]+alpha*h)
f_1 = func(x_1)  
for i in x_1:
    pol_1.append(polynom_value(g, x, i))
for i in range(n+1):
    D_1.append(f_1[i]-pol_1[i])

print("{:^12}{:^12}{:^12}{:^12}".format("x~_i", "f~_i", "pol~_i", "f(x~_i)-pol(x~_i)"))
for i in range (n+1):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_1[i], f_1[i], pol_1[i], D_1[i]))
print()

#all same for Chebyshev nodes
x_2 = []
f_2 = []
pol_2 = []
D_2 =[]

for i in range(n+1):
    x_2.append(0.5*((a+b)+(b-a)*cos((2*i+1)*pi/(2*(n+1)))))
f_2 = func(x_2)
x_2.reverse()
f_2.reverse()
A_2 = matrixA(x_2)
g_2 = Gauss_back(A_2, f_2[:])
g_2.reverse()
x_2.reverse()
f_2.reverse()
for i in x_2:
    pol_2.append(polynom_value(g_2, x_2, i))
for i in range(n+1):
    D_2.append(f_2[i]-pol_2[i])

print("{:^12}{:^12}{:^12}{:^12}".format("x^_i", "f^_i", "pol^_i", "f(x^_i)-pol(x^_i)"))
for i in range (n+1):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_2[i], f_2[i], pol_2[i], D_2[i]))
print()

#comparing polynom's values for charged x and Chebyshev nodes
pol_3 = []
for i in x_1:
    pol_3.append(polynom_value(g_2, x_2, i))

D_3 =[]
for i in range(n+1):
    D_3.append(f_1[i]-pol_3[i])
    
print("{:^12}{:^12}{:^12}{:^12}".format("x~_i", "f~_i", "pol^^_i", "f(x~_i)-pol^^(x~_i)"))
for i in range (n):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_1[i], f_1[i], pol_3[i], D_3[i]))
print()
