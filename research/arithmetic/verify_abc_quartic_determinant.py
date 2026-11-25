from fractions import Fraction as Q
from itertools import product
from math import comb

def determinant(rows):
    m = [row[:] for row in rows]
    sign, prior = 1, 1
    for col in range(len(m)-1):
        pivot = next((i for i in range(col,len(m)) if m[i][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            m[col],m[pivot] = m[pivot],m[col]
            sign *= -1
        val = m[col][col]
        for i in range(col+1,len(m)):
            for j in range(col+1,len(m)):
                num=m[i][j]*val-m[i][col]*m[col][j]
                assert num % prior == 0
                m[i][j]=num//prior
            m[i][col]=0
        prior=val
    return sign*m[-1][-1]

for r,p,inverse in product(range(1,4),(3,5,7),(False,True)):
    mon=[(a,b) for a in range(r+1) for b in range(4*(a if inverse else r-a)+1)]
    q=p**4
    units=[z for z in range(1,q) if z%p][:len(mon)]
    points=[((2*pow(z,-4 if inverse else 4,q))%q,z) for z in units]
    rows=[[x**a*z**b for a,b in mon] for x,z in points]
    det=determinant(rows)
    s=4*comb(r+1,3)+comb(r+1,2)
    assert det != 0 and det % (p**(4*s)) == 0, (r,p,inverse)
print("18 nonzero exact determinants satisfy the prime-power divisibility")

r=100
delta=Q(1,500)
s=4*comb(r+1,3)+comb(r+1,2)
a=r*(r+1)*(8*r+7)//6
b=2*a
margin_inverse=(8-48*delta)*s-a-Q(29,20)*b
margin_positive=(8-48*delta)*s-s-Q(5,3)*b
assert margin_inverse==Q(53833,5)
assert margin_positive==Q(1633574,15)
prime_exponent=(Q(98,25)+2)/12+2*delta
assert prime_exponent==Q(373,750)<Q(1,2)-delta==Q(249,500)
print("Exact margins and both stratum exponents verified")
