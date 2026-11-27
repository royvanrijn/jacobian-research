# Marked three-torsion equations beyond the elliptic sign quotients

The rational points and full Mordell–Weil group are determined in
[HG6](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
This note gives exact polynomial inputs and a local deformation calculation
for the next construction problem. It does **not** supply a geometrically
simple rational example with three independent three-torsion directions,
an infinite family with maximal rational
three-torsion, or a rational-Weierstrass curve with three virtual
directions. HG8 below does give a family transferring unbounded
three-primary Sha, with a genus-two companion whose rank remains uncontrolled.
HG9 explains why the present rational isogeny class cannot realize exact
odd full torsion. HG10 completely classifies one repeated-factor gluing
chart and excludes it from the rational-conic family search.
HG11 supplies a separate infinite family of geometrically simple
Jacobians with exact torsion [3,3] and a rational Weierstrass point.
HG12 imposes simultaneous all-prime virtuality for its two marked
functions. These are intermediate construction results, not solutions
of the three-direction or exact [6,6] targets.
The separate [HG13–HG14 construction](GENUS_THREE_RANK_ZERO_THREE_SHA.md)
now gives rank-zero genus-three Jacobians with arbitrarily large
elementary three-Sha. It uses a different HLP elliptic family; the
specific HG8 companion-rank problem remains open.

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
The [alternative HLP family, HG13–HG14](GENUS_THREE_RANK_ZERO_THREE_SHA.md),
achieves the general rank-controlled Sha existence objective with rank
zero. It does not bound the rank of the specific \(D_R\) above.

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

## A rational two-torsion obstruction to exact odd full torsion

**Proposition HG9.** Let \(A/\mathbb Q\) be isogenous over
\(\mathbb Q\) to a product of elliptic curves each having a nonzero
rational point of order two. Every abelian variety over \(\mathbb Q\)
isogenous to \(A\) has a nonzero rational point of order two.
In particular no member of the rational isogeny class of HG1 has exact
full torsion [3,3,3] or [3,9]. Moreover a rational non-split
(2,2,2) gluing of three elliptic curves cannot yield exact odd full
torsion.

**Proof.** On each elliptic factor the mod-two Galois image fixes a
nonzero vector, so is a subgroup of the order-two stabilizer in
\(\operatorname{GL}_2(\mathbb F_2)\). The kernel of reduction from
\(\operatorname{GL}_2(\mathbb Z_2)\) is pro-two. Consequently the
Galois image on the product's two-adic Tate module is pro-two. A
rational isogeny identifies the rational two-adic representations.
Replacing the integral Tate lattice therefore still gives a finite
two-group as the image on the new abelian variety's two-torsion.
A two-group acting on a nonzero vector space over \(\mathbb F_2\)
has a nonzero fixed vector: its nontrivial orbits on the odd-sized set
of nonzero vectors have even size. This proves the first assertion.

For the last assertion, [HLP, Lemma 13 and the paragraph following its
proof, pp. 20–21](sources/howe-leprevost-poonen-split-jacobians.pdf)
classify the non-split maximal isotropic gluing kernels and show that
descent of such a kernel requires a rational nonzero \(Q_i\in E_i[2]\)
on each factor. Thus the first assertion applies to the glued quotient.
This is a consequence of the classical gluing classification, not a
new obstruction to arbitrary constructions of Jacobians. In particular
it does not exclude the even exact group [6,6], nor varieties outside
this rational isogeny class. The marked-norm deformation problem must
leave that class to reach exact [3,3,3].

## Complete classification of a repeated marked factor

Write the HG1 elliptic family in the coordinate \(U=X+1\) as
\[
 E_h:Y^2=U(U^2+A_hU+B_h),\qquad
 A_h=h^2-3,\quad B_h=2h+3.
\]
For the particular repeated-factor chart set
\[
 h_1=h_2=(b-3)/2,\qquad h_3=(bs^2-3)/2,
 \quad b,s\in\mathbb Q,
\]
using the marked two-torsion point \(U=0\) on all three factors.
Exclude \(b=0,1,9\), \(s=0\), and singular third factors. Put
\(A=(b^2-6b-3)/4\). The nondegenerate branch of HLP equation (10)
with this marking has equation
\[
 2(b^2-6b-3)(b^2s^4-6bs^2-3)
       -s\bigl((b^2-6b-3)^2+64b\bigr)=0. \tag{11}
\]
Indeed with \(d_1=d_2\) the other sign in (10) would give
\(A_3/d_3=1\), hence \(B_3=0\). On the displayed branch,
\(A_3/d_3=(A^2+4b)/(A^2-4b)\). Together with
\(A_3^2-d_3^2=4bs^2\) this gives
\(d_3=s(A^2-4b)/(2A)\) after choosing the sign of \(s\), and
\(2AA_3=s(A^2+4b)\), which is (11).
The division by \(A\) loses no rational point, since its two roots
\(b=3\pm2\sqrt3\) are irrational. Changing the two-torsion markings
is outside this chart.

**Proposition HG10.** The normalization of (11) is the elliptic
curve
\[
             E:v^2=u^3+u^2-24u+36
                    =(u+6)(u-2)(u-3).
\]
It has rank zero and rational group [2,4]. Its only admissible rational
images \((b,s)\) are
\[
                 (6,-3/2),\qquad(-3/5,-3).
\]
For both images the HLP quotient conic has no real point. Thus this
repeated marked factor chart provides no rational-conic gluing and
cannot provide an infinite family of the requested kind.

**Exact proof certificate.** Homogenize (11) with
\((x:y:z)=(1:s:b)\), obtaining the irreducible degree-eight
polynomial \(F_0\) in the retained checker. Apply two quadratic
Cremona substitutions \(q\mapsto M_i(yz,xz,xy)^t\), with
\[
 M_1=\begin{pmatrix}0&0&1\\0&1&1\\1&0&1\end{pmatrix},
 \qquad
 M_2=\begin{pmatrix}-1&-1/2&1/8\\0&1&-3/4\\1&0&1\end{pmatrix}.
\]
After removing exceptional factors \(x^4y^4z^3\) and
\(x^2y^3z^2\), and constant contents, the strict transforms have
degrees five and three. The latter is the smooth cubic
\[
 32x^3+144x^2y+162xy^2-96x^2z-216xyz
                -81y^2z+72yz^2-16z^3=0,
\]
with flex \((0:1:0)\). Conversion with this flex gives the stated
elliptic model. Exact descent proves rank zero; the torsion is [2,4].
The eight rational points map to the following eight points on the
original projective plane model:
\[
\begin{gathered}
 (1:-3:-3/5),\ (0:1:0),\ (1:-3:1),\ (1:1:9),\\
 (1:-3/2:6),\ (1:-1/3:9),\ (0:0:1),\ (1:2:0).
\end{gathered}
\]
The checker extends the composed map at its base points by exact
power-series arithmetic on the smooth cubic. Its three sections have
degree four, so their vanishing orders are at most twelve; precision
twenty suffices to recover their projective limits. Every rational
singular point of the original plane model is outside the admissible
chart, so no admissible rational point is lost on normalization.
Excluding the stated degeneracies leaves exactly the two pairs above.

Finally set \(\Delta=A^2-4b\), \(\Delta_3=A_3^2-4bs^2\),
\(d_3=s\Delta/(2A)\) and \(R=\Delta d_3\). Substitution in
HLP Proposition 14 gives conic coefficients
\(d=e=1/(bs)\), \(f=-1/b\), with the required quartic signs.
Thus the quotient conic is
\[
                     X^2+Y^2-sZ^2=0.
\]
Both admissible values of \(s\) are negative. These conics are
definite and have no real, hence no rational, point. This excludes
their presentation as double covers of a rational parameter line;
it does not deny the existence of the hyperelliptic double covers
of these anisotropic conics.

The [retained exact checker](verify_genus_three_hlp_repeated_factor.py)
passed with SageMath 10.9:
```sh
sage -python research/arithmetic/verify_genus_three_hlp_repeated_factor.py
```
Its scope is this chart, not the full three-factor gluing surface.
The next infinite-family search must use a non-diagonal parameter
curve or a different marking and separately certify rationality of
its conic and nonconstant geometric moduli.

## A rational-Weierstrass family with two marked directions

Two markings can be built without elliptic quotient symmetries. For
\(a\in\mathbb Q\), define
\[
\begin{split}
 q_1&=t^3+at^2-(4a+30)t-1,\\
 q_2&=t^3+(4a+27)t^2+(5a+27)t-4,\\
 H&=\frac73(t^2+t+1),\qquad
 g_1-g_2=H(q_2-q_1),\\
 g_1+g_2&=(q_1^2+q_1q_2+q_2^2)/H.
\end{split}
\]
The last quotient is polynomial. Explicitly,
\[
\begin{split}
g_1={}&\left(\frac72a+\frac{225}7\right)t^4
 +\left(\frac{241}{14}a+\frac{803}7\right)t^3
 +\left(\frac92a^2+\frac{127}2a+\frac{1620}7\right)t^2
 +(6a+45)t+1,\\
g_2={}&\left(-\frac72a-\frac{216}7\right)t^4
 +\left(-\frac{151}{14}a-\frac{569}7\right)t^3
 +\left(\frac92a^2+\frac{71}2a+\frac{297}7\right)t^2
 +(-15a-81)t+8.
\end{split}
\]
Put \(f_a=(g_1^2+q_1^3)/t^2\). It is a monic degree-seven
polynomial in \(\mathbb Q[a,t]\), and
\[
                  g_i^2-f_at^2=-q_i^3\quad(i=1,2). \tag{12}
\]
The factorization of a difference of cubes proves the equality of the
two numerators. Their constant and linear coefficients vanish, proving
divisibility by \(t^2\). This construction keeps the value of each
marked function at \(t=0\) a cube; HG12 uses that feature locally.

**Theorem HG11.** For every integer \(a\equiv-9\pmod{185}\),
the smooth projective curve \(C_a:y^2=f_a(t)\) has genus three,
a rational Weierstrass point at infinity, a geometrically simple
Jacobian, and
\[
                   J_a(\mathbb Q)_{\rm tors}\simeq[3,3].
\]
Infinitely many geometric isomorphism classes occur among these curves.
This theorem supplies two marked directions for the construction
campaign. It does not realize either exact [6,6] or [3,3,3], and no
historical-priority claim is made for this intermediate torsion group.

**Proof.** At \(a=-9\),
\[
\begin{split}
196f_{-9}(t)={}&196t^7-5211t^6+41022t^5+116149t^4\\
              &-69840t^3+13896t^2+4480t-1008. \tag{13}
\end{split}
All coefficient denominators in the family divide powers of 2 and 7.
The specified congruence fixes its reductions at 5 and 37. Both are
smooth, with Jacobian orders 108 and 46215, whose gcd is 9.
The good-reduction torsion bound of [Müller–Reitsma, Section
1.1](sources/muller-reitsma-genus-three-torsion-2211.03372v4.pdf)
therefore bounds the full rational torsion by 9.

Above the zeros of each \(q_i\), select the divisor \(D_i\)
on which \(g_i+ty\) vanishes. Coprimality of \(g_i,q_it\),
verified already in the common reduction at 5, ensures that conjugate
norm factors have no common zero there. Equation (12) gives
\(\operatorname{div}(g_i+ty)=3D_i-9\infty\), since
\(\operatorname{ord}_\infty(t)=-2\),
\(\operatorname{ord}_\infty(y)=-7\), and the \(ty\) term
has pole order nine, greater than the degree-four term's order eight.
Thus \([D_i-3\infty]\in J_a[3](\mathbb Q)\).
The retained checker verifies that their reductions span nine elements
at 5 (and at 37). This proves exact torsion [3,3]. The odd-degree
model gives the rational Weierstrass point.

The Frobenius polynomial at 5 is
\[
                    P(X)=X^6-X^5+8X^3-25X+125.
\]
It is irreducible over \(\mathbb Q\). Moreover no quotient of
two distinct roots is a root of unity. Here is the complete finite
certificate for the latter assertion. Such a quotient has degree at
most 36, so its order \(m\) has \(\varphi(m)\le36\).
The prime-factor formula gives \(\varphi(m)^2\ge m/2\),
hence \(m\le2592\). The polynomial
\[
      \operatorname{Res}_X(P(X),P(zX))/(z-1)^6
\]
has no factor \(\Phi_m(z)\) for any of the 71 integers
\(2\le m\le2592\) with \(\varphi(m)\le36\).
These are exact rational polynomial calculations. Consequently the
six conjugate Frobenius roots remain distinct after every positive
integer power. Each powered polynomial stays irreducible, so the
reduction remains simple over every finite extension. A proper
geometric abelian subvariety in characteristic zero would specialize,
after finite extension, to one in this good reduction. This proves
geometric simplicity for all parameters in the congruence class.

For geometric distinctness, let \(B_a(X,Z)=Z^8f_a(X/Z)\).
Use unnormalized classical transvectants and put
\[
 J_2=(B_a,B_a)_8,\quad G=(B_a,B_a)_4,\quad J_4=(G,G)_8.
\]
The contraction formula shows that these have coefficient degrees
two and four, and determinant weights eight and sixteen. Thus
\(J_4/J_2^2\) is unchanged by a linear change of binary coordinates
and by scaling the binary form, so it is a geometric invariant of
the hyperelliptic curve where defined. The checker obtains
\[
 \left.\frac{J_4}{J_2^2}\right|_{a=-9}
 =\frac{448591166316421718490}{3947463686715586303},\qquad
 \left.\frac{J_4}{J_2^2}\right|_{a=-8}
 =\frac{182782757428770535}{1009280839235802}.
\]
These unequal values prove that the rational function of \(a\)
is nonconstant. Its fibres are finite, so it takes infinitely many
values on the progression, proving the assertion.

## Two simultaneous all-prime virtuality conditions

For (13), write \(F=196f_{-9}\),
\[
\begin{split}
 q_1&=t^3-9t^2+6t-1,&
 G_1&=9t^4-563t^3+342t^2-126t+14,\\
 q_2&=t^3-9t^2-18t-4,&
 G_2&=9t^4+221t^3+1224t^2+756t+112.
\end{split}
\]
Here \(G_i=14g_i\). On \(Y^2=F(t)\), put
\(A_i=196(G_i+tY)\). This is \(14^3(g_i+ty)\), so it
represents the same cubic Kummer class as the original marked function.

**Proposition HG12.** Put
\[
                M=2^7\,3^3\,7^4\,17\,19\,251\,1163.
\]
For every integer \(n\equiv0\pmod M\) for which \(F(n)\)
is not a square, let \(K_n=\mathbb Q(\sqrt{F(n)})\).
Both \(A_i(n)\) are nonzero algebraic integers and have valuation
divisible by three at every finite prime of \(K_n\). At every
prime above \(2,3,7,17,19,251,1163\) they are local cubes.
No independence of their specialized ideal classes is asserted.

**Proof.** The integral identities are
\[
  \operatorname{N}_{K_n/\mathbb Q}A_i(n)=-(196q_i(n))^3,
  \qquad A_1(0)=14^3,\quad A_2(0)=28^3.
\]
The cubics \(q_i\) are irreducible, so neither norm vanishes at
an integer. Direct resultants give
\[
\begin{split}
 \operatorname{Res}(G_1,q_1)&=-3^3\,17^2\,19\,251,\\
 \operatorname{Res}(G_2,q_2)&=-2^8\,3^3\,17\,19\,1163.
\end{split}
\]
Let \(S=\{2,3,7,17,19,251,1163\}\). At a prime above
\(p\notin S\), simultaneous divisibility of \(A_i(n)\)
and its conjugate would force \(p\mid G_i(n),q_i(n)\),
contradicting the resultants. Thus any prime dividing \(A_i(n)\)
outside \(S\) lies above a split rational prime and its conjugate
does not divide \(A_i(n)\). Its valuation equals that of the norm
and is divisible by three.

For \(p\in S\), write \(c_1=14,c_2=28\). Integrality of
\(Y=\sqrt{F(n)}\) and the polynomial coefficients gives
\(A_i(n)-c_i^3\in n\mathcal O_{K_n}\). By the exponents
chosen in \(M\),
\[
 v_p(n)>3v_p(c_i)+2v_p(3).
\]
After multiplying by the ramification index, this is exactly the
Hensel inequality for a root near 1 of
\(X^3-A_i(n)/c_i^3\) in each completion above \(p\).
It proves the local cube assertion and completes the all-prime
valuation argument.

The congruence condition is substantive: at \(n=1\),
\(F(1)=36\cdot2769\), and \(A_2(1)\) has valuations
4 and 5 at the two primes above 2 in \(\mathbb Q(\sqrt{2769})\).
The checker verifies this negative control as well as the norm,
resultant and Hensel bounds. Thus the construction implements
simultaneous local design for two directions. A third independent
norm identity is still required for the original three-direction
class-group objective; virtuality alone also does not establish
independence in a specialized ideal class group.

Both new results have the retained exact replay
```sh
sage -python research/arithmetic/verify_genus_three_two_norm_family.py
```
It passed under SageMath 10.9. The Frobenius computation at 5 is
cross-checked by direct point enumeration over the fields of sizes
5,25,125. The infinite-family, good-reduction specialization and
all-prime arguments are the written proofs above. This remains an
internal execution, without independent replay or external review.

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
HG9 is the written pro-two Galois argument above using the retained
HLP classification. HG10 has the separate exact replay above, including
the proven rank-zero computation and all eight rational images. Neither
result resolves the broad geometrically simple or infinite-family gaps.
HG11–HG12 have the additional replay above. The infinite family there
has rational three-rank two; it does not resolve the maximal rank-three
family gap or the three-direction virtual-unit gap.
