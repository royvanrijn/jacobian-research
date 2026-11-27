# Global octic conductor exclusion

Let \(F=(P,Q)\) be an actual complex polynomial plane Keller map in
a source frame \((x,y)\) with \(P\) monic of positive degree in \(x\).
Its coordinate-order conductor cannot have global \(x\)-degree eight.
In conjunction with the lower-degree exclusions, every noninvertible
map in this frame has even global conductor degree at least ten.
The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the conductor degree on every source line with twice
the defect of its order. The defect is therefore at least five
where the global leading coefficient is nonzero, with only the
retained lower bound one on an arbitrary source line. The
[degree-ten exclusion](GLOBAL_CONDUCTOR_TEN_EXCLUSION.md) raises
the global lower bound to twelve.

The [characteristic-sequence calculation and cubic–quintic exclusion](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md)
and the [all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md)
reduce global conductor degree eight to sequence \((6,4,5)\).
The two degree-two sequences \((9,2)\) and \((9,6,2)\) also have a
[separate period proof](HYPERELLIPTIC_OCTIC_CONDUCTORS_AND_FOUR_PERIOD_RIGIDITY.md).
We exclude it by integrating a polynomial parameter equation in the
actual degree-four and degree-six target generators. The equations
give one rational parameter, whose finite zeros and poles each
force every source collision gap to diverge.

## 1. The reduced generators and their parameter equation

Set \(K=\mathbb C(y)\). Reduce the target generators over \(K\) as
in the characteristic-sequence calculation. Their degrees are four,
six, and their degree-five approximate root is monic of degree
three in the degree-four generator. This is the approximate-root
statement in
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf),
with the hypotheses verified in the characteristic-sequence calculation.

Normalize the generators to be monic in the source variable. The
approximate root has terms only of weights at most twelve, where
the two variable weights are four and six. In fact its monomials
with exponents zero, one and two in the degree-four variable have
weights distinct modulo six, so no higher terms can cancel upon
substitution. Its top terms must be the cube of the degree-four
generator minus the square of the degree-six generator. Below
weight twelve, the permitted monomials are the product of the
generators, the square of the degree-four generator, the generators
themselves and the constant. Completing the square in the degree-six
generator by an affine shear and translating the degree-four
generator eliminate the product and both square terms apart from
the degree-six generator squared. The remaining affine function
of the degree-four generator has source degree at most four.
Consequently the actual transformed pair \(U,V\) satisfies
\(\deg_u(V^2-U^3)=5\). With the polynomial part at source infinity,
this says exactly
\[
 U=u^4+a_0u^2+b_0u+c_0,\qquad V=[U^{3/2}]_{\ge0}.
\tag{1}
\]
Here the depressed quartic frame \(u=x-h(y)\) has \(h\in K\).
The square-root expansion is the monic formal expansion at
\(u=\infty\). The bracket in (1) gives its polynomial part.
The target change has nonzero determinant \(\lambda\in K\).

For the moment write \(U=u^4+a u^2+b u+c\) in this unscaled
frame and put \(\Delta=a^2-4c\). The leading coefficient of
\(W=V^2-U^3\) is \(3b\Delta/8\) at degree five, so \(b\Delta\ne0\).
The parameter compatibility used in the cubic–quintic proof is
\[
 U_u V_y-U_y V_u=\lambda J+k U_uV+\partial_u H(U,V),
 \qquad k=\lambda'/\lambda,\qquad J=\operatorname{Jac}_{x,y}(P,Q).
\tag{2}
\]
Indeed the parameter derivative of the target automorphism, composed
with its inverse, is a target polynomial field of divergence \(k\).
Removing the field \((0,kV)\) makes it Hamiltonian. This proves (2)
over \(K\) and after finite extensions of \(K\).

The left side of (2) has degree at most two in \(u\).
To see this, write \(E=V-U^{3/2}=O(u^{-1})\) at infinity. The
contribution of \(U^{3/2}\) cancels, and \(U_u E_y-U_y E_u\)
has degree at most two. Also \(3U_uV-2UV_u\) has degree at most
two by substituting the same \(E\). The antiderivative \(H(U,V)\)
in (2) has degree at most ten; pole semigroup \(\langle4,5,6\rangle\)
gives a basis \(1,U,W,V,U^2,UW,UV\) for this filtered space.
Its \(UV\) coefficient is forced by the degree-nine term to be
\(-2k/5\). After combining its derivative with \(kU_uV\),
the result has degree at most two. The possible remaining
nonconstant basis terms have derivatives of distinct degrees
\(3,4,5,7,8\), so all their coefficients vanish.

Choose \(\rho\) with \(\rho^{10}=\lambda\). Change the source
coordinate by \(u=\rho z\), and scale the pair by \(\rho^{-4},\rho^{-6}\)
respectively. The target determinant is now one, the original map
in \((z,y)\) has Jacobian \(J\rho\), and the scaled generators
still have form (1). Use \(a,b,c,\Delta,W\) for their scaled
coefficients and expressions. In this frame (2) becomes simply
\[
 U_z V_y-U_y V_z=J\rho .
\tag{3}
\]
Taking its coefficients in degrees two, one and zero gives
\[
 A=b\Delta\in\mathbb C^*,\qquad
 B=\Delta^2-8ab^2\in\mathbb C,\qquad
 G=b^3/4-3Aa/16,\quad G'=J\rho.
\tag{4}
\]
For example the first coefficient is \(3(b\Delta)'/4\);
the degree-one coefficient is \(-3(\Delta^2-8ab^2)'/32\).
The constant coefficient is the displayed derivative of \(G\)
once \(A\) is constant. See the
[exact checker](verify_global_octic_conductor.py).

Although the rescaling was made over \(K(\rho)\), it is rational.
The nonzero invariant satisfies
\[
 b_0(a_0^2-4c_0)=\rho^7 A .
\]
Thus \(\rho^7,\rho^{10}\in K\) and
\(\rho=(\rho^7)^3/(\rho^{10})^2\in K\). In particular \(a,b,c,G\)
are all in \(K\). Equations (4) give the exact one-parameter form
\[
 \Delta=A/b,\qquad
 a=\frac{A^2}{8b^4}-\frac{B}{8b^2},\qquad
 c=\frac{a^2-\Delta}{4}.
\tag{5}
\]

## 2. Collision equations and the finite-fibre obstruction

Let \(z_\pm\) be a distinct matched pair with sum \(s\) and
difference \(g\). Write \(D=z^2+a/2\). The identity
\(U=D^2+bz-\Delta/4\) gives
\[
 V=\tfrac32 U D-\tfrac12 D^3+\tfrac38 b^2.
\]
Indeed this is the polynomial part in (1), by the binomial expansion
of \((D^2+bz-\Delta/4)^{3/2}\) at infinity. Equal \(U\)-values
give \(s(D_++D_-)=-b\), so \(s\ne0\). Equal \(V\)-values
then give \(3U=D_+^2+D_+D_-+D_-^2\). Since
\(D_+-D_-=sg\) and \(D_++D_-=(s^2+g^2)/2+a\), eliminating
between these two equations yields
\[
 s^4+2a s^2-bs+\tfrac32\Delta=0,
 \qquad g^2=\frac{3(\Delta-2bs)}{2s^2}.
\tag{6}
\]
Only necessity for each matching pair is used.

The unit-gap conclusion of the cubic–quintic proof applies in this
class as well. At any finite parameter \(y_0\), at least one generic
matching pair has both original \(x\)-coordinates finite and distinct
in the limit. The exact specialized conductor is nonzero
nonconstant, a generic branch reaches one of its finite roots,
and monicity of the original \(P\) makes its mate finite.
Equal limits would force the vanishing of both original \(x\)
derivatives by their divided differences. Consequently at least
one physical gap \(\rho g\) is a unit at every finite parameter.

Suppose \(b\) has a zero of order \(q>0\) at a finite parameter.
Then the leading terms in (5) are
\(a\sim A^2/(8b^4)\) and \(\Delta=A/b\);
by (4), \(G\sim-3A^3/(128b^4)\), so \(\rho\) has pole order
\(4q+1\). The Newton polygon of (6) gives two sums with
\(b^2s\to \pm i A/2\) and two with \(s^2\sim-6b^3/A\).
For the first two sums the numerator \(\Delta-2bs\) has
leading term \(A(1\mp i)/b\). For the other two, it is dominated
by \(\Delta\). The gap \(g\) has order \(3q/2\) or \(-2q\),
respectively. Every physical gap \(\rho g\) therefore diverges.

If \(b\) has a pole of order \(q>0\), then \(a=O(b^{-2})\),
\(\Delta=A/b\) and \(G\sim b^3/4\), so \(\rho\) has pole
order \(3q+1\). The quartic in (6) has three sums with
\(s^3\sim b\) and one with \(s\sim 3A/(2b^2)\).
For the first three, \(g^2\sim-3s^2\). For the small sum,
\(bs\sim 3\Delta/2\), so \(g^2\sim-4b^3/(3A)\).
Every physical gap again diverges. The valuations in both cases
may be calculated after extending the local parameter field to the
paired roots, so the argument also covers coincident sums and
multiple collision values.

The unit-gap conclusion excludes all these finite zeros and poles.
A nonzero rational \(b\in\mathbb C(y)\) of this kind is constant,
so (4) and (5) make \(G\) constant, contradicting \(G'=J\rho\ne0\).
This closes the \((6,4,5)\) case and therefore the global degree-eight
conductor case.
