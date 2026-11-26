#!/usr/bin/env python3
"""Replay exact exponential margins and weighted-tree controls."""

from fractions import Fraction as F
from math import factorial, prod
def lo(x, degree=55):
    return sum(x**k/factorial(k) for k in range(degree+1))
def hi(x, degree=55):
    return lo(x,degree)+x**(degree+1)/factorial(degree+1)/(1-x/(degree+2))

eb=hi(F(1))
assert eb<F(68,25) and lo(F(1))>F(19,7)
assert max(2*eb/lo(F(27,4)),256*eb/lo(F(27,2)*F(255,256)))<F(1,128)
z=eb/128
assert z/(1-z)<F(37,200) and z/(1-z)**2<F(1,5)
assert lo(F(11))>59001 and F(27)/(lo(F(11))-1)<F(1,2100)
assert 40*11**2/(lo(F(11))-1)<1
assert hi(F(4,2100)+F(1,59000))<F(101,100)
assert eb/48*hi(F(4,2100))<F(1,16)
assert hi(F(1,59000))/16<F(1,15)
t=F(101,600); z=F(1,15)
assert hi(t)-1+z**8/(1-z)<F(37,200)
assert t*hi(t)+z**8*(8-7*z)/(1-z)**2<F(1,5)
assert lo(F(1))*(1-F(37,200))>1+F(37,200)
a=F(59001,59000)/lo(F(27,2))
assert a<F(1,720000) and 784*a<F(108,100000)
assert eb/32*hi(F(297,5000))<F(1,8)
x=F(6,5)/lo(F(11,2))
assert x/(1-x)<F(1,200)
assert hi(F(1,4200)+F(2,59000))<F(101,100)
assert F(68,25)/(1-F(7,76))<3
assert 8*eb*F(297,2500)/(2*(1-F(297,1250)))<F(17,10)
assert 24*eb*F(108,100000)*F(101,100)*hi(F(17,10))/(1-F(297,1250))<F(13,25)
assert lo(F(143,50))>F(125,8) and lo(F(26,25))>F(125,47)
assert 4*F(101,100)*F(6,5)/(1-F(78,125))<13
assert hi(F(2,5))<F(3,2)
assert F(3,20)*hi(1+F(101,1600))<1
assert 2/lo(F(1361,240))+1/lo(F(18))<F(1,130)
assert F(1,129)<F(1,120) and F(2,190)<F(1,55)
assert F(1,60)+F(13,720000)<F(1,55)
assert hi(F(49,720000))<F(21,20)
assert (F(49,720000)+F(2,55))/(1-F(294,720000))<F(1,25)
assert F(14,55)+F(4410,25*720000)<F(1,3)

# Independent series recurrence for weighted ordered tree shapes.
def mul(a,b,n):
    out=[F(0)]*n
    for i in range(min(n,len(a))):
        for j in range(min(n-i,len(b))):
            out[i+j]+=a[i]*b[j]
    return out
for t in (F(1,11),F(1,20),F(1,80)):
    weights=[prod(1+q*t for q in range(u))/F(factorial(u)) for u in range(12)]
    tree=[F(0)]*12
    for k in range(1,12):
        current=[F(1)]
        series=F(0)
        for u in range(k):
            if len(current)>k-1:
                series+=weights[u]*current[k-1]
            current=mul(current,tree,k)
        tree[k]=series
        expected=prod(k+t*l for l in range(k-1))/F(factorial(k))
        assert tree[k]==expected and expected<=3**(k-1)
print('Exact curvature-27 margins and 33 independent weighted-tree checks passed.')
