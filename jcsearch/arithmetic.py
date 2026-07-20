"""CRT, rational reconstruction, and simple multivariate Hensel/Newton lifting."""
import sympy as sp

def crt_pair(a, m, b, n):
    k = ((b-a)*pow(m, -1, n)) % n
    return (a+m*k) % (m*n), m*n

def rational_reconstruct(value, modulus):
    """Return a/b with |a|,|b| <= sqrt(modulus/2), if it exists."""
    bound = int((modulus//2)**0.5)
    r0, r1, s0, s1 = modulus, value % modulus, 0, 1
    while abs(r1) > bound:
        q = r0//r1
        r0, r1, s0, s1 = r1, r0-q*r1, s1, s0-q*s1
    if s1 and abs(s1) <= bound and sp.gcd(r1,s1)==1 and (value*s1-r1) % modulus == 0:
        if s1 < 0: r1, s1 = -r1, -s1
        return sp.Rational(r1,s1)
    return None

def hensel_newton(polynomials, variables, root, prime, levels=3):
    """Lift a nonsingular root mod p to mod p**levels."""
    current=[int(v)%prime for v in root]
    modulus=prime
    jac=sp.Matrix(polynomials).jacobian(variables)
    for _ in range(1,levels):
        values=sp.Matrix([int(f.subs(dict(zip(variables,current)))) for f in polynomials])
        if any(int(v)%modulus for v in values): raise ValueError("not a root at current precision")
        rhs=sp.Matrix([(-int(v)//modulus)%prime for v in values])
        J=jac.subs(dict(zip(variables,current))).applyfunc(lambda e:int(e)%prime)
        delta=J.inv_mod(prime)*rhs % prime
        current=[int(a+modulus*int(d)) for a,d in zip(current,delta)]
        modulus*=prime
    return current, modulus

