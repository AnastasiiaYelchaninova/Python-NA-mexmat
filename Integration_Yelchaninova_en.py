#calculate integral by trapezoid and Simpson methods (using parabolas)
#integral from 0 to 1 of the function (x*ln(x))/(1+x)

from math import *

def bounds():
    a = 0+0.0000000001
    b = 1
    return a,b

def f(x):
    return (x*log(x))/(1+x)

def trapezoid(f,Iold,k):
    a,b = bounds()
    if k == 1:
        Inew = (f(a) + f(b))*(b - a)/2.0
    else:
        n = 2**(k-2) # amount of new points
        h = (b - a)/n # point spacing
        x = a + h/2.0 # first point coordinates
        sum = 0.0
        for i in range(n):
            sum = sum + f(x)
            x = x + h
        Inew = (Iold + h*sum)/2.0
    return Inew

Iold = 0.0
for k in range(1,21):
    Inew = trapezoid(f,Iold,k)
    if (k > 1) and (abs(Inew - Iold)) < 1.0e-6: break
    Iold = Inew

def simpsonsRule(n):
    a,b = bounds()
    sum = float()
    sum += f(a) 
    sum += f(b)
    h=(b-a)/(2*n) # point spacing
    oddSum = float()
    evenSum = float()
    for i in range(1,n): # calculate odd values (except for the first and last)
        oddSum += f(a+2*h*i)
    sum += oddSum * 2
    for i in range(1,n+1): # calculate even values (except for the first and last)
        evenSum += f(a+h*(-1+2*i))
    sum += evenSum * 4
    return sum * h/3
Isimpsons = simpsonsRule(100)

n = 10
Inew1 = simpsonsRule(10)
Iold1 = Inew1 + 1
while abs(Inew1 - Iold1) >= 1.0e-6 :
    Iold1 = Inew1
    n = n*2
    Inew1 = simpsonsRule(n)

Imain = -0.177532966575868
print("Task: to calculate the integral from 0 to 1 of the function (x*ln(x))/(1+x) using the trapezoid and Simpson methods (with parabolas)")
print()
print("Manual integral value = ", Imain)
print("Value of integral by trapezoid method = ", Inew)
print("Simpson integral value = ", Inew1)
print("Calculation errors: ", Imain-Inew, Imain-Isimpsons)
