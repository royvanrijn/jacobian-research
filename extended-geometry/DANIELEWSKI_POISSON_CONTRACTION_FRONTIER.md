# The Danielewski Poisson-contraction obstruction

The multi-dicritical construction leaves one precise affine-target problem.
For a squarefree polynomial

\[
P(x)=xQ(x),\qquad
D_P=\{cw=P(x)\}\subset\mathbb A^3_{c,x,w},
\]

the source map is an open immersion into
\(\mathbb A^1_a\times D_P\).  Replacing the Danielewski factor by two
polynomial target coordinates while preserving \(a\) would turn it into an
ordinary Keller map of affine three-space.

This note reduces that replacement to one Poisson equation.  A period of
the residue form then excludes every polynomial solution when
\(\deg P>1\).  A direct affine-linear calculation supplies an algebraic
cross-check.  Thus a successful affine contraction cannot preserve the
suspension coordinate: it must genuinely mix \(a\) into at least two target
coordinates.

## 1. The intrinsic bracket

On the chart \(c\ne0\), use the residue form

\[
\omega_P=\frac{dc\wedge dx}{c}.
\]

For \(F,G\in B_P=k[c,x,w]/(cw-P(x))\), define the bracket by

\[
\begin{aligned}
\{F,G\}={}&
c(F_cG_x-F_xG_c)\\
&+P'(x)(F_cG_w-F_wG_c)\\
&+w(F_xG_w-F_wG_x).
\end{aligned}                                      \tag{1}
\]

It satisfies

\[
\{c,x\}=c,\qquad
\{c,w\}=P'(x),\qquad
\{x,w\}=w,                                         \tag{2}
\]

and \(cw-P(x)\) is a Casimir, so (1) descends to \(B_P\).  On \(c\ne0\),

\[
dF\wedge dG=\{F,G\}\omega_P.                       \tag{3}
\]

The formula extends across \(c=0\) because (1) is polynomial.

## 2. Exact reduction for an \(a\)-preserving contraction

Recall the open immersion

\[
\Phi_P(a,b,c)=
\left(a,c,bc,bQ(bc)\right)
\colon\mathbb A^3\longrightarrow
\mathbb A^1\times D_P.                             \tag{4}
\]

For \(F,G\in B_P\), form the polynomial map

\[
K_{F,G}=
(a,F\circ\Phi_P,G\circ\Phi_P)
\colon\mathbb A^3\longrightarrow\mathbb A^3.       \tag{5}
\]

Since

\[
\Phi_P^*(da\wedge\omega_P)
=-da\wedge db\wedge dc,
\]

equation (3) gives

\[
K_{F,G}^*(da\wedge dF\wedge dG)
=-\{F,G\}\circ\Phi_P\,
 da\wedge db\wedge dc.                             \tag{6}
\]

Consequently

\[
\boxed{\{F,G\}=1}
                                                        \tag{7}
\]

produces an ordinary polynomial Keller map
\(\mathbb A^3\to\mathbb A^3\) with determinant \(-1\).

This implication does not by itself prove that the missing Danielewski
divisors remain dicritical after (5).  A solution of (7) must still be
checked for nonproperness and for separation of the boundary components.
It is nevertheless the exact algebraic gate for every contraction that
preserves \(a\), or replaces it by \(\lambda a+H(c,x,w)\).

The evident contraction \((c,x)\) does not pass the gate:

\[
\{c,x\}=c.
\]

Its vanishing factor is exactly the factor removed by the residue form.
No polynomial change of coordinates on the target can turn this determinant
into a nonzero constant.

## 3. Pure coordinate obstructions

Write \(\delta_F=\{F,-\}\).

For \(F=c\),

\[
\delta_c(c)=0,\qquad
\delta_c(x)=c,\qquad
\delta_c(w)=P'(x).                                 \tag{8}
\]

This is a locally nilpotent derivation.  After localizing at \(c\), it is
\(c\,\partial_x\), so

\[
\ker\delta_c=k[c].
\]

If \(\{c,G\}=1\), then \(G\) is a slice.  The slice theorem would give

\[
B_P=k[c,G].
\]

But after base change to an algebraic closure, the fiber \(c=0\) of the
left side is

\[
k[x,w]/(P(x)),
\]

which has \(\deg P\) irreducible components, whereas the fiber \(c=0\) of
\(k[c,G]\) is the integral affine line.  This is impossible for
\(\deg P>1\).  The same argument with \(c\) and \(w\) exchanged excludes
\(\{w,G\}=1\).

For \(F=x\), give \(B_P\) the grading

\[
\deg c=-1,\qquad \deg x=0,\qquad \deg w=1.
\]

Then

\[
\delta_x(c)=-c,\qquad
\delta_x(x)=0,\qquad
\delta_x(w)=w.                                     \tag{9}
\]

Thus \(\delta_x\) multiplies each homogeneous component by its weight and
annihilates the weight-zero part.  Its image has zero weight-zero component,
so it cannot contain \(1\).  Hence \(\{x,G\}=1\) is also impossible.

## 4. Direct affine-linear no-go theorem

**Theorem.**  Let \(k\) have characteristic zero, let \(P\) be squarefree
of degree \(n\ge3\), and let

\[
F=\alpha c+\beta x+\gamma w+\delta\in B_P.
\]

There is no \(G\in B_P\) with \(\{F,G\}=1\).

**Proof.**  Constants do not affect \(\delta_F\).  Its values on the three
generators are

\[
\begin{aligned}
\delta_F(c)&=-\beta c-\gamma P'(x),\\
\delta_F(x)&=\alpha c-\gamma w,\\
\delta_F(w)&=\alpha P'(x)+\beta w.                 \tag{10}
\end{aligned}
\]

If \(\beta\ne0\) and \(\alpha\gamma\ne0\), all three expressions vanish at

\[
c=-\frac{\gamma P'(x)}{\beta},\qquad
w=-\frac{\alpha P'(x)}{\beta}
\]

whenever

\[
\alpha\gamma(P'(x))^2-\beta^2P(x)=0.               \tag{11}
\]

The polynomial in (11) has degree \(2n-2>n\), hence is nonconstant and has
a root over an algebraic closure.  At the corresponding point of \(D_P\),
the Hamiltonian vector field \(\delta_F\) vanishes.

If \(\beta\ne0\) and exactly one of \(\alpha,\gamma\) vanishes, choose a
simple root \(r\) of \(P\).  When \(\gamma=0\), the point

\[
(c,x,w)=
\left(0,r,-\frac{\alpha P'(r)}{\beta}\right)
\]

is a zero of (10).  When \(\alpha=0\), use

\[
(c,x,w)=
\left(-\frac{\gamma P'(r)}{\beta},r,0\right).
\]

If both vanish, \(F=\beta x+\delta\), and
\((0,r,0)\) is a zero.

It remains to take \(\beta=0\).  If \(\alpha\gamma\ne0\), choose a root
\(r\) of \(P'\).  Squarefreeness gives \(P(r)\ne0\).  Equations

\[
\alpha c=\gamma w,\qquad cw=P(r)
\]

have a solution over the algebraic closure, and (10) vanishes there.
The cases with only \(\alpha\) or only \(\gamma\) nonzero are the pure
\(c\) and pure \(w\) cases excluded by the slice argument above.

At every critical point \(p\) of \(F\),
\(\{F,G\}(p)=\delta_F(G)(p)=0\) for every \(G\).  It therefore cannot equal
one.  Since a solution over \(k\) would remain a solution after base change,
the result holds over \(k\). \(\square\)

## 5. The residue-period obstruction

The affine-linear theorem is a shadow of a complete obstruction.

**Theorem.**  If \(P\) is squarefree and \(\deg P>1\), there are no
\(F,G\in B_P\) and no \(\lambda\in k^\times\) satisfying

\[
\{F,G\}=\lambda.                                   \tag{12}
\]

**Proof.**  It is enough to work over \(\mathbb C\).  If the data are
defined over an arbitrary characteristic-zero field, they descend to a
finitely generated subfield, which embeds in \(\mathbb C\).

Choose two distinct roots \(\alpha,\beta\) of \(P\), and a simple path
\(x(t)\) between them containing no other root.  Along its interior choose
a continuous square root \(r(t)^2=P(x(t))\), and set

\[
c=r(t)e^{i\theta},\qquad
w=r(t)e^{-i\theta}.
\]

The circle collapses at both endpoints and gives the usual vanishing
two-sphere \(\Sigma_{\alpha\beta}\subset D_P(\mathbb C)\).  On its interior,

\[
\frac{dc}{c}=\frac{dr}{r}+i\,d\theta,
\]

and hence, up to the choice of orientation,

\[
\int_{\Sigma_{\alpha\beta}}\omega_P
=2\pi i(\beta-\alpha)\ne0.                         \tag{13}
\]

Thus \(\omega_P\) is not exact.  But (3) and (12) would give

\[
\omega_P
=\lambda^{-1}dF\wedge dG
=d\!\left(\lambda^{-1}F\,dG\right),
\]

a contradiction. \(\square\)

For \(\deg P=1\), the obstruction disappears exactly as expected:
\(D_P\simeq\mathbb A^2\), and after normalizing \(P'=1\), the pair
\((c,w)\) has bracket one.

## 6. The coupled three-coordinate frontier

For the first arithmetic example

\[
P=x(x^2+1),
\]

the period theorem says that the two geometric dicritical divisors and the
bijective reductions at \(p\equiv3\pmod4\) cannot be transferred to an
affine target by any \(a\)-preserving polynomial contraction, at any
degree.

Write three general target coordinates as

\[
A,F,G\in B_P[a].
\]

Their relative volume coefficient is

\[
\boxed{
A_a\{F,G\}-F_a\{A,G\}+G_a\{A,F\}=1.
}                                                    \tag{14}
\]

Equation (14), not (7), is the remaining affine-target gate.  If at most
one of \(A,F,G\) depends on \(a\), its left side is a product of an
\(a\)-derivative and a bracket of two elements of \(B_P\).  Since the only
units of \(B_P[a]\) are constants, equality to one would force a nonzero
constant bracket, contradicting the period theorem.

Consequently:

\[
\boxed{
\text{any surviving contraction must make at least two target
coordinates depend nontrivially on }a.
}                                                    \tag{15}
\]

For every solution of (14), one must then evaluate pole valuations along
the missing divisors \(c=0,\ x=\alpha\) and retain only those for which the
components remain distinct nonproperness divisors after composition with
(4).

This closes both the full \(a\)-preserving problem and every
one-coordinate mixing.  The unresolved problem is the coupled equation
(14) with at least two active \(a\)-directions.

The first such finite search space is now closed by the
[minimal coupled ansatz theorem](DANIELEWSKI_MINIMAL_COUPLED_ANSATZ.md).
For target coordinates in
\(\langle a,c,x,w,ac,ax,aw\rangle\), every universal Jacobian minor obeys
the coefficient lock
\([b^2c^2]J=3[1]J\); hence no nonzero constant Jacobian is possible, even
when all three coordinates depend on \(a\).

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_poisson_contraction.py
```

The checker verifies the quotient Poisson bracket, both Keller determinant
identities, the Hamiltonian formulas and critical-point constructions, the
locally nilpotent and grading obstructions, the nonzero vanishing-sphere
period, and an exhaustive small-height affine-linear regression for
\(P=x(x^2+1)\).
