# Pure cyclic suspensions of plane polynomials in Keller maps

**Theorem KellerCyclicSuspensionPlaneMate.** Let \(F=(P,Q,R)\) be a
three-variable polynomial Keller map over a characteristic-zero field.
Suppose \(R=Z^m+B(x,y)\) with \(m\ge2\), where \((x,y,Z)\) is a
polynomial source coordinate system of constant nonzero Jacobian.
Then \(B\) has a polynomial Jacobian mate. More precisely, assign
weights \(1,1,k\) to \(x,y,Z\), with \(k\ge1\) and \(\deg B\le mk\).
If the weighted degrees of \(P,Q\) in this frame are at most
\(M_P,M_Q\), then a mate has degree at most \(M_P+M_Q-k\).

This is the cyclic variant of the
[quadratic suspension proof](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md).
We can first extend the field to include a primitive \(m\)-th root
of unity \(\zeta\). On \(Z^m+B=t\), the exact residue identity of
[Wang–Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1)
makes
\[
 \alpha_t=\frac{dx\wedge dy}{mZ^{m-1}}
\]
the derivative of a constant multiple of \(P\,dQ\). The displayed
expression holds on \(Z\ne0\) and is the restriction of the residue
form on the smooth fibre. Pullback by \(Z\mapsto\zeta Z\) multiplies
\(\alpha_t\) by \(\zeta\). Projecting the primitive onto this
character, by the \(\zeta\)-weighted average over the cyclic group,
preserves its derivative.

The \(dx,dy\) coefficients of a character-\(\zeta\) one-form have
\(Z\)-exponents congruent to \(1\) modulo \(m\), and the \(dZ\)
coefficient has exponents congruent to \(0\). Thus reducing powers
by \(Z^m=t-B\) puts the primitive in the form
\(Z\beta_t+D_t\,dZ\). Here \(\beta_t\) is a plane one-form and
\(D_t\) a plane polynomial, both polynomial in \(t\). Subtract
\(d(ZD_t)\) on the fibre, with \(t\) fixed, and replace
\(\beta_t\) by \(\beta_t-dD_t\). Then \(Z\beta_t\) is a primitive.
The relation \(mZ^{m-1}dZ=-dB\) gives
\[
 m(t-B)\,d\beta_t-dB\wedge\beta_t=dx\wedge dy.
\tag{1}
\]
This equality on \(t\ne B(x,y)\) is a polynomial identity in
\(x,y,t\).

If the degree of \(\beta_t\) in \(t\) is \(j>0\), the coefficient
of \(t^{j+1}\) in (1) says that the top coefficient of \(\beta_t\)
is closed on the plane. Write it as \(dH\) by polynomial
integration. The exact form
\[
 d(Z^{m+1}t^{j-1}H)
 = Zt^{j-1}\left((t-B)dH-\frac{m+1}{m}H\,dB\right)
\]
has that same top coefficient after division by \(Z\). Subtract
it to lower the degree in \(t\). If \(\beta_t\) is independent of
\(t\), comparison of the \(t\)-coefficient in (1) makes it closed.
Write \(\beta=dA\). The constant term of (1) is
\(dA\wedge dB=dx\wedge dy\), the claimed Jacobian identity.

Give \(t\) weight \(mk\), and give each differential its coordinate's
weight. Every reduction by \(Z^m=t-B\) preserves the weighted upper
bound. The initial \(\beta_t,D_t\) and \(\beta_t-dD_t\) all have
weight at most \(M_P+M_Q-k\). At each step \(H\) has degree at most
\(M_P+M_Q-k-mkj\), and the subtraction preserves the bound.
The final \(A\) has degree at most \(M_P+M_Q-k\). If an extension
was needed for \(\zeta\), the equation \(\operatorname{Jac}(A,B)=1\)
with this degree bound is a finite system of linear equations over
the original field in the coefficients of \(A\). A solution over
the extension gives a solution over the original field.

## Bounded completions in ordinary degree five

**Corollary QuinticKellerPureCyclicSourceCompletion.** A three-variable
Keller map of ordinary degree at most five over a characteristic-zero
field is invertible if a constant linear combination of its target
components can be written in an affine source frame as
\[
 R=c\bigl(z+s(x,y)\bigr)^m+B(x,y),
 \qquad c\ne0,\quad m\in\{3,4,5\}.
\]
For a quartic component of this form, the possible exponents are
three and four. The inverse-degree bounds for a completion are
six and eight, respectively.

First work over \(\mathbb C\). Put \(Z=c^{1/m}(z+s)\).
For \(s\) nonconstant, the coefficient of \(z\) in \(R\) is
\(mc s^{m-1}\), so
\[
 (m-1)\deg s \le \deg R-1.
\]
The bound also holds for a constant shift. The constant term in
\(z\) is \(cs^m+B\), giving
\(\deg B\le\max(\deg R,m\deg s)\). In ambient degree at most
five, for \(m=3\) this gives \(\deg s\le2,\ \deg B\le6\).
For \(m=4,5\) it gives \(\deg s\le1,\ \deg B\le5\).
If \(R\) is quartic, then \(\deg s\le1,\ \deg B\le4\) for both
admissible exponents.

Take \(k=2\) in the theorem, so \(Z-c^{1/m}s\) has weight at most
two, \(\deg B\le mk\), and \(P,Q\) in the new frame each have
weighted degree at most ten. The plane mate has degree at most
eighteen. The
[retained planar degree bound](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds)
then makes \(B\) a planar coordinate.

Write \(d=\deg B\). If \(d=1\), take an affine mate \(A_0\).
In the inverse of \((A_0,Z,R)\), solve for \(x,y\) from
\((A_0,R-Z^m)\), giving degree at most \(m\). To recover
\(z=Z/c^{1/m}-s(x,y)\) the degree is at most
\(m\max(1,\deg s)\).

If \(d>1\), use
[planar degree reduction](sources/garcia-gwozdziewicz-ploski-coordinate-lines-2015.pdf)
(Proposition 3.4 and Lemma 3.5) as in the quadratic completion.
Take a mate \(A_0\) of degree \(e<d\) with \(e\mid d\), and write
\(B=B_1+H(A_0)\), with \((A_0,B_1)\) a planar automorphism of
degree \(e\) and \(\deg H=d/e\). For \(e=1\), choose a second
affine coordinate in a frame with \(A_0\) first. In this frame
\(B\) is affine linear in the second coordinate with constant
nonzero slope. This gives the same expression with \(B_1\) affine.
The inverse of a planar automorphism has the same degree, as
shown in the
[quadratic completion proof](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
The input to \((A_0,B_1)^{-1}\) is
\((A_0,R-Z^m-H(A_0))\), of degree at most \(\max(m,d/e)\).
Thus the inverse completion has \(x,y\) of degrees at most
\(\max(me,d)\) and \(z\) of degree at most
\(\max(1,\deg s)\max(me,d)\).

For \(m=3\), \(d\le6\) gives \(e\le3\), so the last bound is at
most eighteen. If \(R\) is quartic this becomes at most six
because \(d\le4,\ e\le2\) and the shift is affine. For \(m=4\)
and \(m=5\), the bounds \(d\le5,\ e\le2\) give inverse degrees
at most eight and ten, respectively. The bounds also dominate
the affine-\(B\) case. The
[bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
now applies, since the degree product is at most ninety.
Finitely generated coefficient fields embed in \(\mathbb C\),
and uniqueness of the formal inverse at the origin after
translation descends the resulting polynomial inverse.
