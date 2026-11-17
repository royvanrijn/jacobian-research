"""Exact objects for the cubic adjoint/half-plane separation proof."""
import sympy as s

def tilted(n, gamma, t):
    return s.Poly(sum(s.binomial(n,k)*s.rf(gamma,2*k)*t**k for k in range(n+1)),t)

def gamma_moment(n, alpha, z):
    return s.Poly(sum(s.binomial(n,k)*s.rf(alpha,k)*z**(n-k) for k in range(n+1)),z)

def auxiliary(n, gamma, z):
    coeffs=[sum(s.rf(gamma,k)/s.factorial(k)*s.rf(n+1,j-k)/s.factorial(j-k)/2**(j-k)
                for k in range(j+1)) for j in range(n+1)]
    return s.Poly(sum(s.ff(n,j)*coeffs[j]*z**(n-j) for j in range(n+1)),z)

def reversed_tilted(n, gamma, z):
    return s.Poly(sum(s.binomial(n,k)*s.rf(gamma,2*k)*(-1)**k*z**(2*n-2*k)
                      for k in range(n+1)),z)

def additive_convolution(f,g,z,n):
    return s.Poly(s.expand(sum(s.diff(f,z,k)*s.diff(g,z,n-k).subs(z,0)
                               for k in range(n+1))/s.factorial(n)),z)

def coefficients(n,M,z):
    return (2-2*M)/z,M*(M-1)/z**2-1,2*n/z

def polynomial_operator(B,n,M,z):
    return z*s.diff(B,z,3)-(3*z+M-1)*s.diff(B,z,2)+(2*z+2*M-2)*s.diff(B,z)-2*n*B

def gauged_concomitant(R,B,n,M,z,sign=1):
    # f=exp(-sign*z)*z**(1-M)*B, with its nonzero gauge divided out.
    a,b,_=coefficients(n,M,z);r=-sign+(1-M)/z
    return s.cancel(B*s.diff(R,z,2)+(a*B-s.diff(B,z)-r*B)*s.diff(R,z)
       +(s.diff(B,z,2)+2*r*s.diff(B,z)+(s.diff(r,z)+r*r)*B
           -s.diff(a,z)*B-a*(s.diff(B,z)+r*B)+b*B)*R)
