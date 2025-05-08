from math import*

#note that f(x) = x[i]*log(x[i])+2*x[i]*cos(x[i])

#we'll start from equidistant nodes
a = 1
b = 2
n = 5
x = []
f = []

h = (b-a)/n

def func (x):
    f = x*log(x)+2*x*cos(x)
    return f

for i in range(n+1):
    xi = a + i*h
    x.append(xi)
    fi = func(x[i])
    f.append(fi)

d = len(x)

#calc coefficients of Lagrange polynom

m = []
L = []

def coef(j, f, x): 
    p = 1
    for i in range(d):
        if i == j:
            continue
        else:
            p = p * (x[j]-x[i])
    mi = f/p
    return mi

def Lagr_print(j, x):
    s = ''
    for i in range(d):
        if i == j:
            continue
        else:
            s = s + '(x-' + str(x[i]) + ')*'
    Li = s[: -1]
    return Li

for j in range(d):
    m.append(coef(j, f[j], x))
    L.append(Lagr_print(j, x))

# building Lagrange polynom for equidistant nodes

def Lagrange_polynom(m, L):
    s = 'L(x) = '
    for i in range(len(m)):
        s = s + str(m[i]) + '*' + L[i] + '+'
    return s[: -1].replace('+-', '-')

# function for values of Lagrange polynom at each point

def Lagrange_value(m, x, k):
    s = 0
    for i in range(len(m)):
        p = 1
        for j in range(len(m)):
            if j == i:
                continue
            else:
                p = p * (k - x[j])
        s = s + m[i]*p
    return s

Lval = []
for k in x:
    Li = Lagrange_value(m, x, k)
    Lval.append(Li)

#some extra calc

D =[]
for i in range(n+1):
    Di = f[i]-Lval[i]
    D.append(Di)

# print all data we made: Lagrange polynom, L(x) for x in [a,b]

print("f(x) = x[i]*log(x[i])+2*x[i]*cos(x[i])")
print()
print(Lagrange_polynom(m, L))
print()
print("{:^12}{:^12}{:^12}{:^12}".format("x_i", "f_i", "L_i", "f(x_i)-L(x_i)"))
for i in range (n+1):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x[i], f[i], Lval[i], D[i]))
print()

#recount values of Lagrange polynom for "changed" x

alpha = 0.7
x_1 = []
f_1 = []

for i in range(n+1):
    xi = x[i]+alpha*h
    x_1.append(xi)
    fi = func(x_1[i])
    f_1.append(fi)    

Lval_1 = []
for k in x_1:
    Li = Lagrange_value(m, x, k)
    Lval_1.append(Li)

D_1 =[]
for i in range(n+1):
    Di = f_1[i]-Lval_1[i]
    D_1.append(Di)

print("{:^12}{:^12}{:^12}{:^12}".format("x~_i", "f~_i", "L~_i", "f(x~_i)-L(x~_i)"))
for i in range (n):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_1[i], f_1[i], Lval_1[i], D_1[i]))
print()

#recount values of Lagrange polynom for Chebyshev nodes

x_2 = []
f_2 = []

for i in range(n+1):
    xi = 0.5 *((a+b)+(b-a)*cos((2*i+1)*pi/(2*(n+1))))
    x_2.append(xi)
    fi = func(x_2[i])
    f_2.append(fi)

#we also need to recount all Lagrange polynom for Chebyshev nodes

m_2 = []
L_2 = []
d_2 = len(x_2)

for j in range(d):
    m_2.append(coef(j, f_2[j], x_2))
    L_2.append(Lagr_print(j, x_2))

Lval_2 = []
for k in x_2:
    Li = Lagrange_value(m_2, x_2, k)
    Lval_2.append(Li)

D_2 =[]
for i in range(n+1):
    Di = f_2[i]-Lval_2[i]
    D_2.append(Di)

print("{:^12}{:^12}{:^12}{:^12}".format("x^_i", "f^_i", "L^_i", "f(x^_i)-L(x^_i)"))
for i in range (n+1):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_2[i], f_2[i], Lval_2[i], D_2[i]))
print()

Lval_3 = []
for k in x_1:
    Li = Lagrange_value(m_2, x_2, k)
    Lval_3.append(Li)

D_3 =[]
for i in range(n+1):
    Di = f_1[i]-Lval_3[i]
    D_3.append(Di)

print("{:^12}{:^12}{:^12}{:^12}".format("x^_i", "f^_i", "L^_i", "f(x^_i)-L(x^_i)"))
for i in range (n):
    print("{:^12f}{:^12f}{:^12f}{:^12f}".format(x_1[i], f_1[i], Lval_3[i], D_3[i]))
print()
