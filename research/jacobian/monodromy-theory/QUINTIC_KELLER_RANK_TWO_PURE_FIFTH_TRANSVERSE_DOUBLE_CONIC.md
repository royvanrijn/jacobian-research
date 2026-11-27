# The pure-fifth-power pencil with a transverse double conic

**Theorem QuinticKellerRankTwoPureFifthTransverseDoubleConic.** No
three-variable Keller map of characteristic zero has ordinary
homogeneous leading parts
\[
                  P_5=z^5,\qquad Q_5=z^3q,\qquad R_4=q^2
\]
with \(q\) a geometrically irreducible ternary quadratic whose
restriction to \(z=0\) has rank two.

The proof uses the highest Laurent modes of successive homogeneous
determinant equations. The [pure-fifth-power residual](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_PENCIL_RESIDUAL.md#removing-the-nonzero-coefficient)
is the reason to consider these leading parts.

## The first correction

Embed the finitely generated coefficient field in \(\mathbb C\).
A linear frame preserving \(z\) gives \(q=xy+\gamma z^2\)
with \(\gamma\ne0\). Indeed its rank-two restriction can be made
\(xy\), and adding multiples of \(z\) to \(x,y\) removes the
mixed \(xz,yz\) terms. If \(\gamma=0\), the conic is reducible.

On \(xz\ne0\), use the chart \((z,q,x)\). At \(z=1\), let
\(r=q(x,y,1)\), so \(y=(r-\gamma)/x\).
Write \(H_i\) for the dehomogenization of
\((P_{5-i},Q_{5-i},R_{4-i})\) in \((x,r)\), and set
\[
 H_0=(1,r,r^2),\qquad
 E_i=\operatorname{diag}(5-i,5-i,4-i)H_i,\qquad D=x\partial_x|_r.
\]
The homogeneous determinant of degree \(11-n\), evaluated at
\(z=1\) and up to an overall sign, is
\[
       J_n=\sum_{i+j+k=n}\det(E_i,\partial_r H_j,DH_k).
 \tag{1}
\]
Indeed the determinant of the source change to \((z,q,x)\)
is \(-x\). At \(z=1\) the Euler column \(z\partial_z+
2q\partial_q+x\partial_x\) replaces the \(\partial_z\) column
without changing its determinant, and the factor \(x\) is included
in \(D\). The Keller condition gives \(J_n=0\) for \(n<11\).
The term linear in \(H_n=(p,t,s)\) is
\[
                         D(6r^2p-10rt+5s).
 \tag{2}
\]
This follows from the cross product of the Euler and
\(r\)-derivative columns of \(H_0\). It is also the cofactor row
in the [pure-fifth-power calculation](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_PENCIL_RESIDUAL.md#removing-the-nonzero-coefficient)
at \(\rho=0\).

Here radial means independent of \(x\) in the chart. A Laurent
polynomial is radial exactly when \(D\) kills it. The polynomial
basis \(1,x^i,y^j\) for \(i,j\ge1\) is free over
\(\mathbb C[q,z]\): use \(xy=q-\gamma z^2\) to reduce every
monomial, and use the distinct Laurent modes at \(z=1\) to see
independence on each homogeneous degree.
In \(J_1=0\), reduction of the cofactor row modulo \(r\)
gives \(R_3=qL\) modulo \(\mathbb C[q,z]\), since the coefficient
of \(R_3\) is five and the other entries vanish modulo \(r\).
This also follows by comparing the positive-degree basis coefficients
before taking Laurent modes. A source translation inducing the
linear correction \(-L/2\) to \(q\) removes the nonradial part
of \(R_3\) through \(R_4=q^2\). The quadratic form of \(q\) is
nonsingular, so such a translation exists. It preserves all
leading homogeneous forms.

With \(R_3\) radial, the nonradial coefficients of \(Q_4\) are
\(3r/5\) times the corresponding coefficients of \(P_4\)
by (2). Multiplication by \(r\) costs degree two. Thus the
nonradial part of \(H_1\) is \(U S\), where
\[
              S=(1,3r/5,0),\qquad U=A x^2+B x+C y^2+F y .
 \tag{3}
\]
The coefficients \(A,B,C,F\) are scalars. Here and below the
ordinary degree constraints refer to the original \(x,y,z\).

## Highest modes of the next three equations

The derivations \(D,\partial_r\) preserve Laurent mode at fixed \(r\).
A component of original degree \(d\), dehomogenized at \(z=1\),
has no Laurent mode \(x^\ell\) with \(|\ell|>d\). In \(J_2\),
the mode \(x^4\) cannot come from \(H_2\), whose component
degrees are at most three. It can only come from two copies
of \(A x^2 S\) in \(H_1\), and (1) gives its coefficient
\(-24 A^2r^2/5\). Thus \(A=0\). Swapping \(x,y\) in the
source gives \(C=0\).

With \(U=Bx+Fy\), the known contribution at mode \(x^2\)
of \(J_2\) is \(-12 B^2r^2/5\).
The mode \(x^2\) of \(H_2\) has three scalar coefficients:
\(x^2 r\) would have degree four and cannot occur in a degree-three
component, and \(R_2\) has degree two. By (2), the equation
at this mode gives
\[
                         [x^2]H_2=(B^2/5,0,0).
\]
Consequently the known contribution of mode \(x^3\) to \(J_3\)
comes only from \(H_0\), \(BxS\) in \(H_1\) and
\((B^2x^2/5,0,0)\) in \(H_2\), or from three copies of \(BxS\).
Its coefficient is \(6 B^3r^2/25\). No component of \(H_3\)
can have mode \(x^3\). Hence \(B=0\), and source interchange
gives \(F=0\). We have made \(H_1\) radial.

With \(H_1\) radial, step \(n=2\) gives again the homogeneous
cofactor equation (2). Now \(R_2\) has degree two, so its positive-degree
basis coefficients cannot contain a factor of \(r\), and (2)
first forces it radial. The degree bound on \(Q_3\) then makes
the nonradial part of \(H_2\) equal to \((M x+N y)S\) for scalars
\(M,N\). In \(J_3\), known terms have Laurent modes at most one
in absolute value since \(H_0,H_1\) are radial. Therefore the
mode \(x^2\) of \(H_3\) vanishes by (2) and the independence of
\(1,r,r^2\); the same holds in the \(y\) direction.
In \(J_4\), \(H_4\) has component degrees at most one.
The only contribution to mode \(x^2\) is from two copies of
\(M x S\) at drop two with \(H_0\), giving coefficient
\(-12 M^2 r^2/5\). Thus \(M=0\). Interchanging \(x,y\) gives
\(N=0\), making \(H_2\) radial.

Now apply (2) for \(n=3,4\) with all smaller drops radial.
Reduction modulo \(r\) again makes \(R_{4-n}\) radial.
The degree bound \(\deg Q_{5-n}\le2\) precludes any nonradial
coefficient \(3r/5\) times that of \(P_{5-n}\), so the whole
correction is radial. Constants are already radial. This makes
all three components polynomials in \(q,z\), contradicting the
nonzero determinant.

## Exact mode calculation

The five coefficients used above are short determinant identities.
The following code replays them with exact rational arithmetic.
Each list contains only the pieces that can contribute to the
stated Laurent mode.

~~~python
from itertools import product
from sympy import Matrix, expand, symbols

x, r, A, B, M = symbols("x r A B M")
top = Matrix([1, r, r**2])
zero = Matrix.zeros(3, 1)
direction = Matrix([1, 3*r/5, 0])

def determinant_at_drop(pieces, drop):
    weights = (5, 5, 4)
    euler = [
        Matrix([(weight-i)*value for weight, value in zip(weights, piece)])
        for i, piece in enumerate(pieces)
    ]
    radial = [piece.diff(r) for piece in pieces]
    modes = [x*piece.diff(x) for piece in pieces]
    return expand(sum(
        Matrix.hstack(euler[i], radial[j], modes[k]).det()
        for i, j, k in product(range(len(pieces)), repeat=3)
        if i+j+k == drop
    ))

assert determinant_at_drop([top, A*x**2*direction], 2) == -24*A**2*r**2*x**4/5
assert determinant_at_drop([top, B*x*direction], 2) == -12*B**2*r**2*x**2/5
correction = Matrix([B**2*x**2/5, 0, 0])
assert (determinant_at_drop([top, B*x*direction, correction], 3)
        == 6*B**3*r**2*x**3/25)
assert determinant_at_drop([top, zero, M*x*direction], 4) == -12*M**2*r**2*x**2/5
cofactor = Matrix([5, 5*r, 4*r**2]).cross(top.diff(r))
assert cofactor == Matrix([6*r**2, -10*r, 5])
~~~

## Exclusion of the double-conic boundary

**Corollary QuinticKellerRankTwoDoubleConicBoundaryExclusion.** The
unique cancelling quartic of a noninvertible three-variable
quintic Keller map cannot have geometrically irreducible double-conic
leading form. In particular it has a repeated line.

By the [quintic rank restriction](QUINTIC_KELLER_RANK_ONE_FIXED_DOUBLE_LINE_NORMAL_FORM.md#consequence-for-quintic-inversion)
the top rank is two and the cancelling quartic has nonreduced
leading form. In the double-conic case the
[pure-fifth-power reduction](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_PENCIL_RESIDUAL.md#removing-the-nonzero-coefficient)
gives \(P_5=z^5,Q_5=z^3q,R_4=q^2\).
The conic has a nonzero restriction to \(z=0\), of rank one or two.
The tangent rank-one restriction is excluded by the
[tangent determinant calculation](QUINTIC_KELLER_RANK_TWO_TANGENT_PURE_FIFTH_DOUBLE_CONIC.md)
and the rank-two restriction by the theorem above. A nonreduced
ternary quartic with no repeated line can only be an irreducible
double conic. This gives the claimed repeated-line condition.
