# Marked three-torsion equations beyond the elliptic sign quotients

The rational points and full Mordell–Weil group are determined in
[HG6](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
This note gives exact polynomial inputs and a local deformation calculation
for the next construction problem. It does **not** supply a geometrically
simple rational example, an infinite family with maximal rational
three-torsion, or a rational-Weierstrass curve with three virtual
directions. HG8 below does give a family transferring unbounded
three-primary Sha, with a genus-two companion whose rank remains uncontrolled.

## Primary comparison and the actual target

[Müller–Reitsma, Theorem 1.1](sources/muller-reitsma-genus-three-torsion-2211.03372v4.pdf)
lists the possible exceptions [3,3,3], [3,9], [2,4,4], [6,6] to its
existence statement for **exact full torsion** on geometrically simple
genus-three hyperelliptic Jacobians. Remark 1.2 also distinguishes split
examples. These are the exceptions in that inspected 2023 statement;
the bounded literature check here is not a claim that all four remain
unresolved in later literature. A simple example merely containing
\((\mathbb Z/3)^3\) would not settle the exact [3,3,3] case unless its
other torsion were excluded. HG6 has full torsion [3,6,12], so it does
not realize that exception even without the simplicity condition.

[Howe–Leprévost–Poonen, Proposition 14 and Section 4.2](sources/howe-leprevost-poonen-split-jacobians.pdf)
require chosen rational two-torsion on three elliptic curves, a square
discriminant product and vanishing twisting factor. In their notation
\(E_i:y^2=x(x^2+A_i x+B_i)\), \(\Delta_i=A_i^2-4B_i\),
\(d_i^2=\Delta_i\), \(R=d_1d_2d_3\), the last condition is
\[
 R\left(A_1^2/\Delta_1+A_2^2/\Delta_2+A_3^2/\Delta_3-1\right)
       -2A_1A_2A_3=0.
\]
Their equation (10), \(A_3/d_3=(A_1A_2\pm4\sqrt{B_1B_2})/(d_1d_2)\),
is an equivalent useful gluing condition. Their Proposition 14 presents
a double cover of a conic; rationality of that conic must be imposed
separately when an equation over \(\mathbb Q(t)\) is desired. Varying
inside this construction retains a split Jacobian. It therefore addresses
the infinite split-family question, not geometric simplicity.

[Lupoian, arXiv:2210.02225v4, Section 2](sources/lupoian-genus-three-three-torsion-2210.02225v4.pdf)
uses polynomial identities \(g^2-f\ell^2=cq^3\) to describe three-torsion.
The following lower-degree identities are derived directly on our curve;
no completeness assertion about that paper's normalized torsion charts is
needed. They avoid putting any sign involution into the defining equations.

## Three independent norm identities

Use \(F\) from HG1 and the marked rational point \(\infty_+\). For
\(a=-1,0,2\), put \(q_a=(t-a)^2\) and use this table:

| \(a\) | \(\ell_a\) | \(c_a\) |
|---|---|---:|
| -1 | \(t^2+5t+7\) | -300 |
| 0 | \(t^2+2\) | -200 |
| 2 | \(t^2-4t-2\) | 24 |

\[
\begin{split}
g_{-1}&=t^6+2t^5-13t^4-52t^3-49t^2+38t+22,\\
g_0&=t^6-3t^5-3t^4-12t^3+6t^2-12t-8,\\
g_2&=t^6-7t^5+5t^4+20t^3+50t^2-52t+40.
\end{split}
\]
All three satisfy exactly
\[
                 g_a^2-F\ell_a^2=c_aq_a^3.
\]
The functions \(g_a+\ell_aR\) have divisors
\(6P_a-6\infty_+\), where \(P_{-1}=(-1,9)\),
\(P_0=(0,4)\), \(P_2=(2,36)\). Thus
\(2[P_a-\infty_+]\) are three-torsion classes. Reduction at 7
certifies that their span has 27 elements, so they form a basis of the
same full three-torsion subgroup as the HG1 pullbacks. This changes the
equations used to mark the torsion, not the rational three-torsion space.

Allow \(f,\ell_i,g_i,q_i\) to be monic of degrees 8,2,6,2 respectively,
and let each \(c_i\ne0\) vary. Do **not** retain the double-root
condition on \(q_i\). Require
\[
                      g_i^2-f\ell_i^2-c_iq_i^3=0 \quad(i=1,2,3).
\]
There are 41 coefficients: eight of \(f\), and eleven per row
(two of \(\ell_i\), six of \(g_i\), two of \(q_i\), one of \(c_i\)).
The leading degree-12 coefficients cancel identically, leaving 36
coefficient equations, of degrees zero through eleven in \(t\).

On the open set where \(f\) is squarefree and
\(\gcd(g_i,\ell_iq_i)=1\), these identities yield rational
degree-two effective divisors \(D_i\) with
\(\operatorname{div}(g_i+\ell_iR)=3D_i-6\infty_+\).
Indeed the two conjugate factors have no common zero above a root of
\(q_i\), so their multiplicities are multiples of three. The norm has
degree six; at \(\infty_+\) the pole has order six, and the norm's
nonzero leading coefficient shows no zero or pole at \(\infty_-\).
Consequently \([D_i-2\infty_+]\in J[3](\mathbb Q)\).

## Exact tangent calculation HG7

**Proposition HG7.** At the three displayed rows, the coefficient
Jacobian has size \(36\times41\) and rank 36 over \(\mathbb Q\).
Its kernel has dimension five, and projection onto the eight coefficients
of \(f\) has rank five. Thus this marked-norm scheme is smooth of
dimension five at the fixed example. This dimension counts coefficients;
it is not a dimension claim for an unmarked moduli space.

For a row, the linearized equation is
\[
 2g_i\,\delta g_i-\ell_i^2\delta f
 -2f\ell_i\delta\ell_i-q_i^3\delta c_i
 -3c_iq_i^2\delta q_i=0.
\]
These expressions define the exact matrix without symbolic differentiation
or numerical tolerances. The
[retained checker](../../unpublished/genus-three-rational-three-torsion/verify_arithmetic.py)
also tests persistence of each known base involution,
\[
 \frac{2t+2}{t-2},\qquad\frac{2-t}{t+1},\qquad-\frac2t.
\]
It permits each involution to move as \((at+b)/(ct-a)\) and permits its
multiplier on the binary octic to vary. The linearized condition is
\(\delta(F\circ M-\lambda F)=0\). In every case, only a
three-dimensional subspace of the five-dimensional marked-norm tangent
space satisfies this condition. Since a finite union of proper rational
linear subspaces cannot cover that space, rational tangent directions
exist that break all three known involutions to first order.

Smoothness supplies formal arcs over \(\mathbb Q[[s]]\) through these
directions. On a formal neighbourhood the three classes remain
independent: three-torsion is finite étale in characteristic zero and its
independence locus is open and closed. Neither a formal arc nor an étale
coefficient chart produces a nonconstant rational curve or infinitely
many rational specializations. Losing these involutions also does not
exclude higher-degree elliptic factors or a split abelian-surface factor.

The next precise construction task is to find an algebraic rational
curve on these equations, or a positive-rank genus-one curve with a
nonconstant map to them, then exhibit a rational smooth specialization
and prove its geometric simplicity, for example by an absolutely simple
good reduction. For the exact Müller–Reitsma [3,3,3] target, additionally
prove the full rational torsion is exactly 27. The full-rank tangent
calculation removes a local coefficient obstruction; these global tasks
remain open.

## Infinite split families, virtual units and the real obstruction

For the separate split-family route, impose HLP's discriminant-product
and twisting-factor equations together with a rational point on the
quotient conic. A rational point on that parameter variety alone is
insufficient for an infinite-family theorem. One needs a rational or
positive-rank genus-one parameter curve and a nonconstant moduli
invariant; repeatedly parametrizing the fixed conic does not change the
geometric isomorphism class of the genus-three curve.

HG3 proves that on the original curve's quadratic specializations only
one marked direction is virtual at every prime. The proposed new family
must impose simultaneous valuation divisibility at every bad prime;
renaming the three generators cannot increase the virtual subspace.
For class-group counting one must also prove independence after
specialization and account for the unit kernel and the counting theorem's
hypotheses. Three rational Jacobian torsion classes alone do not give
three independent ideal classes.

There is a further constraint on the suggested rational Weierstrass
point. HG2 proves \(F>0\) on the entire real projective parameter line.
The homogenized binary octic is strictly positive on the real unit
circle, a compact set, so all sufficiently small real coefficient
deformations remain positive. They have no real Weierstrass point and
hence no rational one. A smooth real path cannot change the number of
real branch points. The rational-Weierstrass design target therefore
requires a different real component, or passage through a singular fibre;
it cannot occur in the immediate real neighbourhood just calculated.
This does not obstruct rational points of the curve away from its branch
locus, nor does it obstruct the geometrically simple torsion target.

## A compatible Mordell cover and the remaining rank problem

[MCT15–MCT20](MORDELL_RAY_CASSELS_TATE_PAIRING.md) concern
\(E_R:y^2=x^3-2R^2\). HG5 transfers three-primary arithmetic for an
existing prime-to-three isogeny; it does not construct that isogeny for
new elliptic factors. In particular the fixed HG1 factors have nonzero
\(j\)-invariants, unchanged by quadratic twists, whereas this Mordell
family has \(j=0\). The current three-parameter \(h_i\) curves also
carry rational two-torsion, at \((-1,0)\). For the squarefree products
of inert primes used in MCT, \(2R^2\) is not a rational cube (an odd
prime divisor has valuation two), so \(E_R\) has no rational two-torsion.
It cannot simply replace an HG1 factor in HLP Proposition 14 with its
chosen-rational-two-torsion hypotheses. Its rational three-torsion is
also absent; keeping the three HG1 markings is a separate restriction
unnecessary for a large-Sha theorem.

An elliptic factor together with a genus-two companion avoids that
chosen-two-torsion constraint. The geometry is classical:
[Shaska, arXiv:1306.5284v2, Theorem 1(ii)–(iii) and proof](sources/shaska-genus-three-split-jacobians-1306.5284v2.pdf)
describes the two quotients of an even octic using \((x^2,y)\) and
\((x^2,xy)\). Its normalized classification is over an algebraic
closure; the following direct formulas establish the needed descent
over \(\mathbb Q\) without making a normalization assumption there.

**Theorem HG8.** For every \(R\in\mathbb Q^\times\), let
\[
 \begin{split}
 C_R &: y^2=(x^2-1)\bigl(1-2R^2(x^2-1)^3\bigr),\\
 D_R &: z^2=u(u-1)\bigl(1-2R^2(u-1)^3\bigr),\\
 E_R &: Y^2=X^3-2R^2.
 \end{split}
\]
Then \(C_R,D_R\) are smooth hyperelliptic curves of genera three and
two. The first has rational Weierstrass points \((\pm1,0)\). There are
degree-two maps
\[
 \alpha:C_R\longrightarrow E_R,\qquad
 (X,Y)=\left(\frac1{x^2-1},\frac{y}{(x^2-1)^2}\right),
 \qquad
 \beta:C_R\longrightarrow D_R,\quad(u,z)=(x^2,xy).
\]
Pullback and pushforward give dual isogenies of degree eight between
\(J_R=\operatorname{Jac}(C_R)\) and
\(E_R\times\operatorname{Jac}(D_R)\), whose composites are multiplication
by two. Hence, for every \(n\ge1\),
\[
 \begin{split}
 \operatorname{Sel}_{3^n}(J_R)&\simeq
 \operatorname{Sel}_{3^n}(E_R)\oplus
 \operatorname{Sel}_{3^n}(\operatorname{Jac}(D_R)),\\
 \Sha(J_R)[3^n]&\simeq\Sha(E_R)[3^n]\oplus
                         \Sha(\operatorname{Jac}(D_R))[3^n].
 \end{split}
\]
In particular, using MCT16, for every \(r\ge1\) there are infinitely
many eligible parameters \(R\) for which
\((\mathbb Z/3)^{2r}\) is a subgroup of \(\Sha(J_R)[3]\).
No uniform bound on \(\operatorname{rank}J_R(\mathbb Q)\) is asserted.

**Proof.** With \(a=2R^2\), the discriminants of the displayed degree-eight
and degree-five polynomials are respectively
\[
                  186624a^5(a+1),\qquad -27a^2(a+1)^2.
\]
They are nonzero for every nonzero rational \(R\). Direct substitution
proves both maps. For either map the target abscissa has degree two on
the source's rational \(x\)-line, so comparing the two hyperelliptic
double covers shows that the map has degree two. The sign involution
\(\sigma:(x,y)\mapsto(-x,y)\) fixes \(\alpha\) and acts on \(\beta\)
as the hyperelliptic involution of \(D_R\). It therefore acts by opposite
signs on the two pullback images. The cross push–pull homomorphisms
vanish in the torsion-free Hom group, while diagonal push–pull is [2].
Thus \(\Psi\Phi=[2]\); finite kernel, equality of dimensions and
surjectivity give \(\Phi\Psi=[2]\). Duality gives
\(\deg\Phi=\deg\Psi=2^3=8\).

The HG5 Kummer argument applies with 2 in place of 8 and an elliptic
times genus-two product in place of three elliptic factors. Multiplication
by two is invertible on every three-power torsion group and modulo
\(3^n\) on local point quotients, at every place. This proves both
isomorphisms, and likewise the full three-primary Sha isomorphism,
without a finiteness assumption. MCT16 supplies infinitely many \(R\)
with \(\dim_{\mathbb F_3}\Sha(E_R)[3]\ge2r\); the direct summand gives
the claimed subgroup. This implication does not use the additional
converse theorem in MCT20.

The rank formula makes the remaining task precise:
\[
 \operatorname{rank}J_R(\mathbb Q)=\operatorname{rank}E_R(\mathbb Q)
                +\operatorname{rank}\operatorname{Jac}(D_R)(\mathbb Q).
\]
MCT16 bounds the first term by one on its nonsingular targets. The
second term must now be bounded on an infinite such subfamily to obtain
the requested controlled-rank theorem. The construction is split and
does not preserve three independent rational three-torsion classes.
The unrestricted geometric simplicity target remains separate.

The [exact symbolic checker](verify_genus_three_mordell_sha_cover.py) verifies
the discriminants with indeterminate \(a\), both maps, their degrees and
signs, and the singular negative control at \(a=-1\). Its command is
```sh
sage -python research/arithmetic/verify_genus_three_mordell_sha_cover.py
```
The run passed under SageMath 10.9. It does not replay MCT16 or compute
the companion rank, Selmer groups or Sha groups. The result is an explicit
classical quotient construction applied to the existing Mordell theorem,
not a new general gluing method or a priority claim for large Sha.

## Verification scope

The replay command is the HG6 command, SageMath 10.9 with exact rational
and finite-field arithmetic and no random input. Besides the arithmetic
theorem, it checks the three norm identities, their reduction independence,
the coefficient-Jacobian rank and all three moving-involution tangent
dimensions. The formal-deformation interpretation and the distinction
between formal and rational families are written arguments. Independent
replay, external mathematical review and formal verification are absent.
HG8 has its separate symbolic replay above; its unboundedness conclusion
retains MCT16's recorded proof and review scope.
