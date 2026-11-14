# A radial boundary budget for the singleton-zero Keller gate

13 November 2026. HC4MYPP1 is a written geometric deduction; HC4MYPP2
is an elementary control family with a written all-genus argument.
Neither result uses Tang's structural theorem. Work over C for topology;
the blowup and valuation arguments hold over any algebraically closed
field of characteristic zero. This is an application of standard
normalization, valuation and pencil techniques, without a priority claim.

## 1. The hypotheses and what a section actually proves

Let G=(G1,G2):A2 -> A2 satisfy

\[
 (G_1,G_2)=(x,y),\qquad J(G)=c\in\mathbb C^*.
\]

Let d=[C(x,y):C(G1,G2)]. The retained relative-linear row orbit has c=-2;
the arguments here do not depend on that normalization. The ideal
condition gives one reduced source point over zero. The previous
[normalization gate](../hc4-relative-row-orbit/PROOF.md) proves that a
nonautomorphic such G has d>1 and omitted normalization length d-1 at zero.

Put X=Bl_0(A2_source), Y=Bl_0(A2_target). Since G is etale, it is flat.
Its pullback of the target maximal ideal is exactly the source maximal
ideal. Flat base change for blowups therefore gives a cartesian square

\[
 X=\mathbb A^2_{\rm source}\mathbin{\times}_{\mathbb A^2_{\rm target}}Y
 \longrightarrow Y.
\]

The lifted map Gtilde is etale. It maps the exceptional E_source=P1
isomorphically onto E_target=P1 by the projectivized linear map DG(0).
In particular it is an isomorphism on the completed local rings along
that source section, but this says nothing about omitted branches.
The base-change input is [Stacks tag0805](https://stacks.math.columbia.edu/tag/0805).

For a direct chart check, write (x,y)=(a,a*b),
G1=a*A(a,b), G2=a*B(a,b), and use the target chart (G1,G2/G1).
The lifted determinant is c/A. Near any point of E_source where that
chart applies, A is a unit. The other target chart treats the remaining
directions. The checker verifies the general chain-rule identity and a
nonlinear automorphism example, without treating this as global
properness evidence.

The projection Y -> P1 is the smooth projection of the total space of
O_P1(-1). Hence pi:X -> P1 is a smooth surjective morphism with the
exceptional curve as a section. Its fibre over a line ell through zero
is isomorphic to G^-1(ell): blowing up its single smooth point at zero
does not alter the curve. Every member is therefore a smooth affine
curve, though a special member may have several disjoint components.

The generic member is geometrically irreducible. Here is the global step
behind that assertion. Resolve the pencil on a smooth projective rational
surface S containing X, and take the Stein factorization of S -> P1.
The exceptional section gives a section of its finite curve factor
T -> P1. For integral T, degrees of that section and the finite map
multiply to one, so the finite map has degree one. The generic projective
fibre is geometrically connected, and is smooth in characteristic zero;
thus it is geometrically irreducible. Removing boundary points preserves
that property. This proves generic irreducibility, not rationality or
irreducibility of every special member.

## 2. All omitted radial divisors must be counted

Let Z -> Y be the finite normalization in C(x,y). Z is a normal surface;
it is finite flat of rank d over smooth Y, by the same normal-surface
flatness argument used in HC4MYOR2. Zariski Main gives an open immersion
X -> Z compatible with the maps to Y. The exceptional source section is
proper over C, so its image in separated Z is closed as well as lying
in that open set. It is an entire prime divisor D0 above E_target.

Write the divisorial pullback as

\[
 E_{\rm target}|_Z=D_0+\sum_{i=1}^k e_iD_i.
\]

For each other prime divisor D_i above E_target, let

\[
 f_i=[\mathbb C(D_i):\mathbb C(E_{\rm target})].
\]

The source section has e0=f0=1. Every other D_i is omitted from X:
the full preimage of E_target inside X is E_source, and an open subset
meeting D_i would meet its generic point. At the generic point of the
target exceptional curve, the integral closure of its discrete valuation
ring gives the standard degree identity d=sum(e_i*f_i). Consequently

\[
 \boxed{\sum_{i=1}^k e_i f_i=d-1.}\tag{R}
\]

This is [the DVR degree formula, Stacks tag0EXQ, Remark15.113.6](https://stacks.math.columbia.edu/tag/0EXQ),
with the actual source contribution removed. It refines the previous
missing-length count into radial divisors and their degrees. In
particular d>1 forces at least one omitted divisor dominating P1.
An incomplete list of ramified divisors cannot certify (R): terms with
e_i=1 are permitted and must also be counted. No assertion that every
omitted divisor is ramified is made.

Choose a general pencil direction. Avoid the finite sets containing
singular points of Z, branch values of the finite D_i -> E_target and
intersections of distinct D_i. The normalized curve over that direction
then has, over radial value zero, the source point of multiplicity one
and f_i omitted points of multiplicity e_i for each i. In this precise
sense (R) is a boundary budget in every general radial slice. It is
stronger data than just the presence of an exceptional section.

The remaining geometric imports, already pinned for HC4MYOR2, are
[Zariski Main](https://stacks.math.columbia.edu/tag/02LR) and
[finite-to-smooth surface flatness](https://stacks.math.columbia.edu/tag/00R3).
The absence of every omitted D_i implies d=1, hence invertibility of
the Keller map. Neither the section nor the ideal hypothesis proves
that absence. We have not classified Z for an unrestricted candidate.

## 3. An independent Euler budget on a general pencil curve

Let S_G be the nonproper-value curve if G is nonproper, and write its
irreducible components as S_j. Define delta_j=d-n_j, where n_j is the
number of source points over a general point of S_j. Let m_j be the
projective degree of S_j and mu_j its multiplicity at zero, taking
mu_j=0 when zero is not on S_j. These are complete curve and source-fibre
data, not merely a discriminant multiplicity.
That S_G is an algebraic curve, with a covering away from it in the
etale case, is the standard
[Jelonek nonproperness theorem](https://www.impan.pl/en/publishing-house/journals-and-series/annales-polonici-mathematici/all/58/3/107876/the-set-of-points-at-which-a-polynomial-map-is-not-proper),
Ann. Polon. Math.58(1993),259-266. The finite generic strata of each
component follow from constructibility of a finite-type morphism.

Let C_ell=G^-1(ell) for a sufficiently general line through zero. All
nonzero intersections with S_G are then distinct generic smooth points
of its components; their number on S_j is m_j-mu_j. Avoiding finitely
many further directions handles all exceptional finite strata. Over
ell minus S_G the restriction is finite etale of degree d. Additivity
of compactly supported Euler characteristic, or the covering-space
cell decomposition, gives

\[
 \chi(C_\ell)=d-\sum_{a\in\ell\cap S_G}
       (d-\#G^{-1}(a)).
\]

If d>1, the zero term is exactly d-1. Therefore, if g is the genus of
the smooth projective completion of C_ell and r its number of punctures,

\[
 \boxed{\sum_j\delta_j(m_j-\mu_j)=2g+r-1.}\tag{E}
\]

For d=1 the map is an automorphism, S_G is empty and C_ell=A1, so the
same final equality reads 0=0. The generic irreducibility proved above
is needed to use this single g and r rather than a componentwise sum.

If d>1, the omitted radial points from (R) and the poles of the radial
function are distinct punctures. Thus

\[
 r\ge1+\sum_i f_i\ge2,\qquad
 \sum_j\delta_j(m_j-\mu_j)\ge2g+\sum_i f_i\ge1.
\]

In particular an alleged nonautomorphic pointed Keller candidate with
zero nonzero-radial deficit budget is impossible. A proposed boundary
ledger must satisfy both (R) and (E); passing them is only necessary.
The genus and puncture count have not been determined for general G.

## 4. What the published pencil theorem would add

[Nguyen Van Chau, Ann. Polon. Math.101(2011),47-53](https://www.impan.pl/shop/en/publication/transaction/download/product/85398)
assumes **all** curves aG1+bG2=0 are irreducible and rational. Its main
theorem then gives invertibility when zero is a regular value, in
particular for a Keller map. Its pencil lemma, in our one-simple-basepoint
case and assuming a rational generic member, reads
sum_lambda(r_lambda-1)=h_infinity-1. Neither rationality nor the vanishing
of that reducibility sum follows from our row-orbit theorem. We do not
import the stronger assertion suggested by omitting those hypotheses.
The final journal theorem is Theorem1; the arXiv v3 labels it Theorem2.

## 5. Same ideal and same Jacobian, every generic pencil genus

Here is an explicit calibration entirely inside HC4MYOR1's oriented
orbit. For each integer g>=0 set

\[
 H_g=(y+xy^2+x^{2g+1},\ 2x),\quad
 D_g=\begin{pmatrix}1&0\\(y^2+x^{2g})/2&1\end{pmatrix}.
\]

Then det(D_g)=1 and (y,2x)*D_g=H_g. Thus H_g=gamma*C0*D_g,
with C0 the retained normalizer. Explicitly

\[
 x=H_{g,2}/2,\quad
 y=H_{g,1}-(y^2+x^{2g})H_{g,2}/2,
\]

so every row has zero ideal (x,y). Every one also has the **same**
Jacobian polynomial

\[
 J(H_g)=-2(1+2xy),\qquad J(H_g)(0)=-2.
\]

These rows are not Keller maps: the common critical curve is xy=-1/2.
Keeping that limitation is essential.

On the pencil chart H_{g,1}+t*H_{g,2}=0, put W=2xy+1. Multiplication
by 4x gives the exact identity

\[
 W^2=1-4x^{2g+2}-8tx^2.\tag{H}
\]

It is a birational equation for the member, with inverse y=(W-1)/(2x)
on x!=0. The right side has degree 2g+2 and is generically squarefree:
at t=0 it is 1-4x^{2g+2}, which is squarefree in characteristic zero.
Consequently (H) is a geometrically integral double cover of P1_x with
2g+2 simple finite branch points and no branch at infinity. The
Riemann-Hurwitz formula gives genus exactly g, for every g>=0.

For completeness, the affine member has precisely three punctures on
this general smooth completion. At x=0, W=1 gives the original origin:
W-1 is divisible by x^2, so y is regular there. The other point W=-1
has a pole in y and is omitted. The two points above x=infinity are
also omitted, since x is an affine coordinate. There are no other
punctures, as the inverse is regular when x!=0. Thus r=3 and b1=2g+2.

The single family keeps the oriented ideal, full Jacobian polynomial
and critical curve fixed while changing generic pencil genus without
bound. Even an elementary relative matrix therefore does not preserve
the rationality needed by the published pencil theorem. The Keller
condition might impose further restrictions; this control does not
settle what those are. It prevents a proof from obtaining rationality
from ideal equality or determinant-one row operations alone.

## 6. A check on the retained nontrivial residual row

There is a similar concrete control already present before any new
matrix choice. For gamma=(-yP(v),xQ(v)), v=xy, the pencil
gamma1+t*gamma2=0 has the birational presentation

\[
 (yP(v))^2=t\,vP(v)Q(v).
\]

The checker proves v*P*Q squarefree of degree eight, hence this general
member has genus three. Yet gamma*C0=(y,2x) has genus-zero pencil.
Also J(gamma)=(v*P*Q)' evaluated at v=xy, so gamma itself is not Keller.
This is a second control using precisely the inherited row; it does
not change the scope of the elementary family in Section5.

## 7. A complete radial slice of the known source-A3 control

The Euler calculation itself applies to any finite-degree map when its
line preimage is a curve. It can therefore be calibrated on the retained
[F1 marked-root example](../hc4-relative-row-orbit/PROOF.md), although
that three-variable map is not a plane candidate. Take the target line
(a,b,c)=(t,t,t). In its finite cubic normalization the root equation is

\[
 t(\zeta^3+\zeta-2)-2\zeta^2=0,\qquad
 t=\frac{2\zeta^2}{D},\quad D=\zeta^3+\zeta-2.
\]

Thus the completed line curve is P1_zeta, a cover of degree three. Put
q=zeta^3-zeta+4. The derivative is -2*zeta*q/D^2. Its four simple
critical points are zeta=0 and the three roots of q. The roots of D
are three simple poles, disjoint from those four points. The source
retains precisely simple roots in the cubic normalization, so its
curve is P1 minus those seven points. Infinity is retained: it gives
the source origin. Hence g=0 and r=7.

For a direct reconstruction, source coordinates along that line are

\[
 X=\frac D{\zeta q},\quad
 Y=\frac{2\zeta(\zeta-3)}D,\quad
 Z=\frac{2X-3X^2Y-t}{X^3}.
\]

Literal substitution gives F1(X,Y,Z)=(t,t,t), and Y+1/X=zeta.
All three coordinates tend to zero as zeta tends to infinity. X has
poles at zero and the three roots of q; Y has poles at all three
roots of D. The complete simple-projective-root correspondence is the
previously proved F1 result, not inferred merely from this substitution.

The three nonzero critical values are the roots of

\[
 R(t)=28t^3-18t^2-t+16,
\]

since Res_zeta(q,tD-2*zeta^2)=-8R(t). They are distinct and nonzero.
At each of these values the cubic has one double and one simple root;
the source fibre has one point, so delta=2. Zero has the same deficit:
its double root is omitted and its simple projective root is the origin.
Consequently the total deficit is 2+3*2=8 and

\[
 \chi(C)=3-8=-5=2-7,\qquad
 \sum_{a\ne0}\delta(a)=6=2g+r-1.
\]

This is an exact nonproper calibration of (E), using a fully described
curve rather than chosen abstract deficit numbers. It does not apply
the plane-specific surface proof of (R) to a threefold. The known
source-A3 counterexample and its full cubic model retain their prior
attribution; the radial slice is a follow-up calculation.

## 8. Verification boundary and next mathematical obligation

The checker replays literal polynomial identities, the lift determinant,
the retained degree-eight squarefreeness, and representatives g=0..8 of
the hyperelliptic control. Arbitrary genus follows from the written
squarefreeness and Riemann-Hurwitz argument, not that finite list.
Normalization existence, the open immersion, the DVR degree formula,
generic irreducibility and Euler additivity are written geometric steps.
Their source pins and the exact finite replay have separate roles.

The next useful input is a **complete** finite-normalization or resolved
pencil model of an actual Keller row beyond the known output-degree
gate. It must identify every omitted radial divisor, including e=1
ones, and supply the curve deficits, genus and punctures for (E).
Neither budget constructs such a row or excludes all of them. The
singleton-zero problem, general mixed/coisotropic route and JC2 remain
open. No formal verification, external review or novelty priority is
claimed.
