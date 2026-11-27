# Three conductor-eighteen exclusions by parameter equations

## The sequences \((14,4,7)\) and \((15,6,4)\)

**Theorem ConductorEighteenSquareRootLaurentExclusions.** Let
\(F=(P,Q)\) be a complex polynomial plane Keller map with \(P\)
monic in the longitudinal source coordinate \(x\). Its reduced
generic directional order cannot have either of the characteristic
sequences \((14,4,7)\) and \((15,6,4)\).

Put \(K=\mathbb C(y)\). Write \(n,m,r\) for the three entries
of the sequence, and \(V,W\) for the reduced generators of
source degrees \(n,m\). By
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf),
the approximate root of source degree \(r\) is monic of
degree \(e=2n/m\) in \(W\). The values of \(e\) in the two
cases are seven and five. The weights \(mi\) for
\(0\le i<e\) are distinct modulo \(n\). Hence terms above
weight \(2n\) cannot cancel in the approximate root. With both
generators source monic its weight-\(2n\) terms are \(W^e-V^2\).
The remaining terms involving \(V\) are linear in \(V\) with
coefficients polynomial in \(W\). Completing the square in \(V\)
and translating \(W\) to depress the right side give
\[
 V^2=A(W)+\eta T,\qquad
 A(t)=t^e+\sum_{j=0}^{e-2}c_jt^j,\qquad
 \eta\ne0,\quad \deg_x T=r .
 \tag{1}
\]
These target changes are triangular or affine and preserve
the stated generator degrees.

Let \(\lambda\in K^\times\) be the determinant of the resulting
target change from \(P,Q\) to the ordered pair \((W,V)\).
Set \(x=\rho z\) with \(\rho^{m+n}=\lambda\) and multiply
\(W,V\) by \(\rho^{-m},\rho^{-n}\), respectively. The monic
generators in \(z\) have the form (1) with rescaled coefficients
and target determinant one. The
[parameter identity](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch)
gives
\[
 W_zV_y-W_yV_z=J\rho+\partial_z H(W,V),
 \tag{2}
\]
where \(J\ne0\) is the original constant Jacobian and \(H\)
is a polynomial over the extended parameter field.

In the monic Laurent branch at infinity, \(A(W)^{1/2}\)
has leading term \(z^n\). Equation (1) gives
\[
 E=V-A(W)^{1/2}
   =\frac{\eta T}{V+A(W)^{1/2}}=O(z^{r-n}).
\]
The degrees of \(W_z,W_y\) are at most \(m-1\).
Thus \(W_zE_y-W_yE_z\) has degree at most \(m-1+r-n\).
This is \(-4\) in the first sequence and \(-6\) in the
second, so it has zero polynomial part. Here degree of a
Laurent series means the highest exponent in \(z\).

The remaining Laurent Jacobian follows by the chain rule:
the contributions from differentiating \(W\) inside \(A(W)\)
cancel, and the original polynomial Jacobian is
\[
 W_zV_y-W_yV_z=
 \left[\frac{W_z}{2 A(W)^{1/2}}\sum_{j=0}^{e-2}c_j' W^j
 \right]_{\ge0}.
 \tag{3}
\]
The term indexed by \(j\), if \(c_j'\ne0\), has Laurent
degree \(mj+m-1-n\) and leading coefficient \(m c_j'/2\).
All terms of negative degree contribute zero. The remaining
terms and the successive semigroup comparisons are as follows.
At each row of the table the higher-index terms have already
been killed. The degree bound on \(H(W,V)\) follows from
(2). Its actual upper degree is the largest semigroup value
within that bound.

| Sequence | \(j\) | Leading degree in (3) | Bound on \(\deg_z H(W,V)\) | Largest allowed degree |
|---|---:|---:|---:|---:|
| \((14,4,7)\) | 5 | 9 | 10 | 8 |
| \((14,4,7)\) | 4 | 5 | 6 | 4 |
| \((14,4,7)\) | 3 | 1 | 2 | 0 |
| \((15,6,4)\) | 3 | 8 | 9 | 8 |
| \((15,6,4)\) | 2 | 2 | 3 | 0 |

In the first case the semigroup is \(\langle4,7\rangle\);
in the second it is \(\langle4,6,15\rangle\). Each row
forces \(c_j'=0\) because the degree of \(\partial_z H\)
is strictly below the indicated leading degree. Once
these terms vanish, all the terms of (3) have zero
polynomial part. Equation (2) then forces \(H(W,V)\) to
have degree one, which neither semigroup admits.

## The sequence \((8,6,7)\)

**Theorem ConductorEighteenCubicApproximateRootExclusion.**
Under the same hypotheses the reduced generic directional
order cannot have characteristic sequence \((8,6,7)\).

Let the reduced generators \(W,V\) have degrees six and eight.
The degree-seven approximate root is monic of degree four
in \(W\). The weights \(6i\) for \(0\le i<4\) are distinct
modulo eight, so terms of weight above 24 cannot cancel.
After making the generators source monic, the top terms
in the approximate root are \(W^4-V^3\). The coefficient
of \(V^2\) is at most affine in \(W\) and that of \(V\)
has degree at most two in \(W\). Depress the cubic in \(V\)
by a triangular shear in \(W\). Then translate \(W\) to
depress the polynomial on the other side:
\[
 V^3+A(W)V=B(W)+T,\qquad
 A(t)=a_2t^2+a_1t+a_0,\qquad
 B(t)=t^4+b_2t^2+b_1t+b_0,\qquad \deg_x T=7 .
 \tag{4}
\]
The degree-seven \(T\) is the approximate root up to a
nonzero scalar. These changes preserve monicity and the
degrees of both generators.

Set \(x=\rho z\), where \(\rho^{14}\) is the determinant
of the target change from \(P,Q\) to \((W,V)\). Rescale
\(W,V\) by \(\rho^{-6},\rho^{-8}\). Equation (4)
retains its form over the parameter extension, and the
target determinant is one. Thus (2) applies to this pair.

Use the derivation \(W_z\partial_y-W_y\partial_z\) on
(4). It annihilates \(W\) and differentiates each target
coefficient by its \(y\)-derivative times \(W_z\), giving
\[
 (3V^2+A(W))\operatorname{Jac}_{z,y}(W,V)
   =W_z\bigl(B_y^{\mathrm{coeff}}(W)
             -A_y^{\mathrm{coeff}}(W)V\bigr)
      +\operatorname{Jac}_{z,y}(W,T) .
 \tag{5}
\]
The superscript indicates differentiation of coefficients
while holding the argument fixed. The denominator
\(3V^2+A(W)\) has leading term \(3z^{16}\). The last
term on the right has degree at most twelve:
\(W_z,W_y\) have degrees at most five, and
\(T_y,T_z\) have degrees at most seven and six. Thus it
has no contribution to the polynomial part after division
by the denominator.

The coefficient-derivative terms on the right of (5)
have distinct degrees whenever they are nonzero.
After division, the only terms whose Laurent degrees
can be nonnegative have leading degrees \(9,3,1\),
from \(a_2',a_1',b_2'\), respectively. Their leading
coefficients are \(-2a_2',-2a_1',2b_2'\) by monicity
and \(W_z=6z^5+\cdots\). For example \(W_zW^2V\)
has degree \(5+12+8=25\), giving the first
degree-nine term. The other coefficient derivatives
give degrees \(-3,-5,-11\) after division.

The semigroup of the reduced order is
\(\langle6,7,8\rangle\). If the Jacobian in (5) has
degree at most nine, (2) bounds \(H(W,V)\) by degree
ten. The largest nongap through ten is eight, so
\(\partial_z H\) has degree at most seven, forcing
\(a_2'=0\). Now the Jacobian has degree at most three,
and the bound \(\deg_z H\le4\) forces \(H\) to be
constant in \(z\). The coefficients in degrees three
and one give \(a_1'=b_2'=0\). All remaining terms in
(5) give zero polynomial part after division. Equation
(2) would then give \(J\rho=0\), a contradiction.
