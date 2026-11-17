"""Van Hoeij--Vidunas F11 over Q(zeta_10); exact polynomial data."""
import sympy as s
a,x=s.symbols("a x")
minimal=a**4-a**3+a**2-a+1
alpha=s.CRootOf(minimal,0)
K=s.QQ.algebraic_field(alpha)
def poly(expr):
    # Reduce in the formal coefficient variable before selecting an embedding.
    expanded=s.Poly(s.expand(expr),x)
    coefficients={m:K([s.QQ.from_sympy(v) for v in s.Poly(s.rem(c,minimal,a),a).all_coeffs()])
                  for m,c in expanded.terms()}
    return s.Poly.from_dict(coefficients,x,domain=K)
scale=s.Rational(64,605)*(-21115+22305*a-51661*a**2+15789*a**3)
A=poly(3025*x**2+(157+2276*a-4828*a**3-2799*a**2)*x-3506*a-1362-2632*a**3+1049*a**2)
B=poly(x**2-(a**2+1)*x+a+a**3)
C=poly(44*x**2+(22*a-88*a**3-44*a**2)*x-42-54*a+9*a**2-16*a**3)
E=poly(2*x-2*a**3-a)
N=poly(scale)*A*B**5
D=C**5*E
# Source F11 is phi=N/D. Move its order-two branch value to infinity.
# beta=phi/(phi-1)=N/(N-D).
Q=N-D
