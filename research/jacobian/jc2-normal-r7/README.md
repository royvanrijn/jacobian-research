# The normal r=7 affine cusp gate: two classes and a conditional exclusion

**Complete coefficient family:** [PF2R7FAM1/PF2R7FAM2](../jc2-r7-image-family/README.md)
now describe every normalized gap-21 image by two explicit rational
one-parameter charts, with all nonimmersion parameters identified. This
applies with additional ramified primes allowed. It supplies the first
image equations as a family; the second image and global matching remain open.

**Affine strata complete:** [PF2R7KCOLL1/PF2R7KCOLL2](../jc2-r7-laurent-collisions/README.md)
complete the Laurent collision analysis alongside polynomial L. Both generic
opens have twelve nodes; all268 nonnodal coefficient points and their affine
types are explicit. [PF2R7TOP1](../jc2-r7-complement/README.md) now
proves proper marked comparison on both nodal opens and at most92
individual-first-image complement types. Compatible additional images
and their joint configuration remain required. The same proof directory
gives a conditional degree(3,5) second-image contact bound and pencil.

**Current single-prime closure:** [PF2OMIT5](../jc2-omitted-divisors/R7_INFINITY_COMPLETION.md)
excludes the entire declared degree-six, single ramified e2,f1 gap21 stratum.
The [reciprocal-contact certificates](../jc2-r7-reciprocal/README.md), combined
with the universal quartic result, give an additional algebraic route.
Any surviving degree-six candidate containing this image requires a second
ramified target image. Its global equations and attachments remain open;
[PF2OMIT6/PF2OMIT7](../jc2-omitted-divisors/README.md) already fix the immersed
residual residue degrees, local lengths and singleton radial packets.

**Current two-class result:** [PF2R7EX2](../jc2-r7-quartic-contact/README.md)
now excludes both nonimmersed classes under the explicit entire single
ramified e2,f1 prime hypothesis **without source finality**. It uses
degree-six conic coefficients, a cubic exclusion and exact quartic-contact
certificates. Earlier finality-dependent and local statements below are
preserved proof history. The immersed and other ramification patterns remain.

13 November 2026. Registered result: **PF2R7CUSP1**. This is an exact
reduction of one of the surviving gates in the
[degree-six F2 Stein programme](../../foundations/plane-jc/F2_GEOMETRIC_DEGREE_SIX_STEIN_REDUCTION.md).
The normal r=7 row, the normal r=5 row, the nonnormal conductor families,
the pair (75,125), and JC(2) remain unresolved.

**14 November scoped closure:** [PF2R7EX1](../r7-adjoint-exclusion/README.md)
excludes both nonimmersed coefficient classes below under PF2R7BND1's
explicit unique ramified affine-prime hypothesis. Its source-model
Borisov inputs and written geometric reduction are stated separately from
the exact algebraic certificates. The immersed family and other possible
ramification patterns remain unresolved. The historical local and
attachment limitations below describe the earlier stages of this route.

**Whole-family extension:** [PF2R7PAR2/PF2R7INF1](../jc2-r7-universal-contact/README.md)
now cover the immersed gap-21 family too. Without source finality they
exclude balanced infinity indices (3,3) and leave three explicit horizontal
infinity budgets. With the additional named Borisov finality input they
exclude the entire single-ramified-prime family. The full normal-r7 row
still requires an audit of that entire ramification hypothesis.

The subsequent [global boundary audit](../jc2-boundary-attachment/README.md),
PF2R7BND1, strengthens the local excess test using Borisov's plane theorem.
Under the single ramified-row hypothesis, either nonimmersed class requires
two or three additional unramified boundary primes with immersed images
different from its cusp curve. The locally allowed A5 collision is globally
incompatible with the boundary. Extra nonproper images and the immersed
r7 family remain unexcluded; this is not a closure of normal r7.

The [minimal-resolution positivity follow-up](../minimal-resolution-positivity/R7_WINDOW.md)
now forces general-line genus two, Picard rank 8--11, four possible generic
infinity decompositions and two to seven contracted boundary primes under
that same single-ramified-row hypothesis. It retains the additional
unramified images and does not close either coefficient class.

## The result

Let p,q be a birational polynomial parametrization of degrees (6,10) in
characteristic zero, with first odd infinity gap 21. Up to the normalization
operations below, every affine nonimmersion is represented by exactly one
of the following two coefficient rows:

| row | a | b | c | d | B | D | E | complete affine singularities |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| I | -4/3 | -2/3 | 4/9 | 4/9 | 7/9 | 496/2187 | 2072/19683 | A2 + A5 + 8 A1 |
| II | -9/5 | -4/15 | 629/675 | 92/675 | -7/9 | -1927/6075 | 41636/273375 | A2 + 11 A1 |

The parametrizations are defined without omitted coefficients by

\[
p=t^6+a t^4+b t^3+c t^2+d t,
\qquad
q=[p^{5/3}+B p^{4/3}+D p^{2/3}+E p^{1/3}]_+.
\]

The bracket means polynomial part at infinity, expanded with leading
coefficient one. Full expanded rational equations are retained in
[`classification.json`](classification.json). Each pair is birational, has
gap exactly 21, and has precisely one nonimmersion parameter, t=1, which is
an ordinary A2 cusp. These are coefficient-normalization classes, not a
claim of inequivalence under every polynomial automorphism of the plane.

**Conditional local exclusion:** if the complete affine dicritical ledger
has just one (e,f)=(2,1) row, row II is impossible by the local excess
identity. Completeness of that missing-boundary ledger is not established
by this follow-up. Row I can meet that numerical identity and is not
excluded here. In addition, no gap-21 pair has a fiber
containing three distinct smooth branches with the same vertical tangent.

## Normalization and complete critical-point classification

Use affine translation of the parameter to center p, translation/scaling
of p and q, and the target shear q -> q+C p+F to reach the displayed form.
These operations preserve nonimmersion. The
[gap theorem PF2D6O1](../../foundations/plane-jc/F2_DEGREE_6_10_POLYNOMIAL_GAP.md)
gives exact odd-tail obstructions rho1,rho3,rho5,rho7,rho9,rho11. Here

\[
\rho_1=\rho_3=\rho_5=\rho_7=\rho_9=0,
\qquad \rho_{11}\ne0.
\]

The recurrence generating these polynomials is reused as a dependency,
not reconstructed from a missing old artifact. The literal identity

\[
(\rho_5+2a\rho_3/3+c\rho_1/3)|_{b=0}=-5d^3/81
\]

shows that b=0 would force d=0 and both p and q even. Birationality excludes
that case, so b is nonzero.

First exclude a critical point at the centered parameter value zero. Scale
the parameter so b=1. The equations p'(0)=q'(0)=0 give d=0 and

\[
D=2Ba^2/9-2Bc/3-10a^3/81+5ac/9+5/54.
\]

The equation rho1=0 then solves for E with constant nonzero coefficient.
After these substitutions let f3,f5,f7,f9 be the positive-denominator-cleared
remaining obstructions, exactly as in
[`critical-zero-elimination.json`](critical-zero-elimination.json). They obey
the compact identity

\[
1=\frac{96a^4-120ca^2+36c^2-52a}{7}f_3
 +\frac{-40a^3+20ca-3}{21}f_5
 +\frac{-232a^2+54c}{63}f_7-\frac{10a}{63}f_9.
\]

Thus every critical parameter is nonzero. Scale it to 1, preserving the
monic/centered normalization. The equation p'(1)=0 eliminates
d=-6-4a-3b-2c; q'(1)=0 eliminates E with coefficient 2. Introduce z=1/b.
The remaining five obstructions and bz-1 define an ideal in
Q[D,B,c,a,b,z]. It is exactly the triangular ideal

\[
\begin{aligned}
0={}&(2z+3)^3(4z+15),\\
0={}&405b+32z^3+264z^2+756z+918,\\
0={}&2799360a+216272z^3+1234044z^2+1153656z+3416283,\\
0={}&125971200c-33224464z^3-224297388z^2-421411032z-295567191,\\
0={}&3359232B+733712z^3+6275052z^2+14525784z+7533351,\\
0={}&5101833600D+1361840944z^3+9638006148z^2+19490644872z+10989597861.
\end{aligned}
\]

The scheme has length four and support z=-3/2,-15/4. Its nonreduced
multiplicity-three point is retained; it is not treated as three additional
geometric classes. Substitution gives rows I and II.

Both ideal inclusions are certified. Forward inclusion is literal
substitution into Q[z]/((2z+3)^3(4z+15)). Reverse inclusion is a rational
six-by-six polynomial lift matrix generated by Singular, with all six
columns independently replayed by FLINT multiplication and addition.
No Gröbner membership decision is used in that replay.

The default Singular lift was unexpectedly large: 496,667,014 bytes. It is
preserved in the ignored local evidence cache, not placed in the source tree:

`artifacts/local/jc2-normal-r7-20261113/r7-critical-one-lift.log`

SHA-256:
`334b52837abb2c679cfe2ce1dba6290d00d34c37e068f1165c35dd159cc6a731`.
[`lift_receipt.json`](lift_receipt.json) and the six `lift-checked-N.json`
receipts bind it to the small generator and target files. Each independent
column replay takes approximately two seconds in this environment. Missing
local evidence does not authorize automatic regeneration. The retained
`r7-critical-one-lift.sing` is the explicit regeneration input; the generation
finished under its 30-second cap, but a future run should enforce an output
size cap or use a more economical representation.

## Complete affine singularity calculation

For each displayed pair, the exact normalization gcd is u-t and
gcd(p',q')=t-1. At t=1, p'' is nonzero and
p''q'''-q''p''' is nonzero, proving an ordinary A2 cusp. Its image fiber
gcd is exactly (t-1)^2, so no additional normalization branch meets it.

Form the exact conductor resultant

\[
R(t)=\operatorname{Res}_u
\left(\frac{p(u)-p(t)}{u-t},\frac{q(u)-q(t)}{u-t}\right).
\]

Up to a nonzero scalar its factorizations are

\[
R_{\rm I}=(t-1)^2(3t^2+3t+1)^3R_{16},
\qquad R_{\rm II}=(t-1)^2R_{22}.
\]

R16 and R22 are squarefree and coprime to the displayed special factors.
For row I the quadratic roots are smooth normalization points sharing
exactly the image (-1/9,49/59049); the image-fiber gcd is the monic quadratic
itself. Conductor order three on each branch gives mutual intersection
three, hence A5. The remaining simple conductor roots pair into ordinary
nodes. Thus the complete affine types are A2+A5+8A1 and A2+11A1. Their
delta totals are both 12, agreeing with the infinity calculation.

## The conditional excess test and the missing-boundary caveat

The canonical Stein proof, equation (3.49), gives

\[
5=\sum_l e_lf_l+\sum_{l,x}\epsilon_{l,x},
\qquad\epsilon_{l,x}=\mu_x-e_lq_x\ge0.
\]

Assume that the complete affine dicritical ledger consists of the simply
ramified component with e=2,f=1, so it requires total excess 3. Its normalization map has residue
degree one, so each affine normalization point has q_x=1.

At a node, local meridian transpositions commute, so every touched local
orbit has size two and contributes zero. At an ordinary cusp, the two
fiber meridians satisfy aba=bab. Their transposition images are equal or
intersecting; the touched orbit has size two or three. The cusp is one
normalization branch, so it contributes at most 3-2=1.

Under that complete-ledger assumption, row II has total excess at most one,
contradicting three. Row I
has a two-branch A5 singularity as well. Its transpositions may generate
S3, contributing two units, so 1+2=3 is possible. Both singularities would
have to attain their maxima in any surviving cover. The finite permutation
checks are retained in [`local-excess.json`](local-excess.json).

It would be incorrect to infer that assumption merely from a classification
of ramification. For the normal cubic map

\[
(\pi,w)\longmapsto(\pi,z=w^3+\pi^7w),
\]

the discriminant pullback factors exactly as

\[
-4\pi^{21}-27z^2
=-(3w^2+\pi^7)^2(3w^2+4\pi^7).
\]

The first curve is ramified; the second is an unramified companion. Deleting
only the first curve, or deleting both, gives two etale dense opens with the
same finite normalization, normal terminal slice and ramification pattern,
but different missing-boundary ledgers. Removing pi=0 in both opens gives
the same example over the affine side of the terminal boundary. Neither
open is asserted to be A2. Thus this is a local insufficiency example, not a
counterexample to a global Keller theorem. A global source-A2 argument is
still required to determine whether the companion is missing. This
follow-up does not certify the completeness step in the older degree-six
ledger or invalidate all of its results.

## Excluding a vertical common-tangent triple

Three distinct smooth branches in one p-fiber with vertical tangents would
force p-P0 to have three double roots. Since deg p=6, it is the square of
a monic cubic. With centered p and b=1 this forces

\[
p=(t^3+ut+1/2)^2-1/4,
\quad (a,c,d)=(2u,u^2,u).
\]

For all three branches to have one target image, the coefficients of t and
t^2 in q modulo that cubic must vanish. Together with the five odd-tail
equations, the seven denominator-cleared polynomials f0,...,f6 in
[`vertical-triple-input.json`](vertical-triple-input.json) satisfy

\[
1=(4u^2/7)f_0+(4/21)f_2
  -((4u^3+1)/21)f_5-(4u^2/21)f_6.
\]

This is checked by direct multiplication. It excludes that projection-specific
triple-tangent gate, not all possible immersed curves or all tangent slopes.

## Bounded topology work and the remaining gate

For the old sharp immersed witness p=t6+t3,
q=t10+5t7/3+7t4/9+7t/81, Sage/SIROCCO produced six certified braid words.
The [complete polynomial-family collision classification](../jc2-r7-collisions/README.md)
now identifies this as the special d=0 member with three ordinary triple
points and three nodes. Its saved case label `generic` is historical;
it is not generic in the complete coefficient family.
An exhaustive S6 transposition enumeration, modulo simultaneous conjugacy
of the first meridian, leaves exactly one assignment and no transitive one.
That excludes this particular witness, not an entire coefficient component.

The cusp calculations did not finish within the stated 60/120-second caps;
no cusp braid certificate or cusp monodromy exclusion is claimed. The
finite-slope three-branch common-tangent ideal was frozen in
`finite-triple-tangent.sing`; a 30-second slimgb attempt returned no completed
basis. `halt 1` is an interruption, not a unit certificate.

The most concrete next tasks are auditing the complete missing-boundary
ledger, a complete local-excess/immersed-family reduction, a smaller exact finite-slope common-tangent elimination, and a
conditioned or partially certified braid analysis of row I. None is replaced
by a generic-witness inference.

## Reproduction and assurance

```sh
.venv/bin/python research/jacobian/jc2-normal-r7/verify_r7.py
.venv/bin/python research/jacobian/jc2-normal-r7/local_excess.py
.venv/bin/python research/jacobian/jc2-normal-r7/check_monodromy.py
```

The first command replays the large local certificate in six checkpointed
workers. `--compact-only` checks small identities and previously hash-bound
receipts but explicitly does not claim a new replay of that certificate.
The algebra, normalization gcds, jets, resultants, and local finite
permutation tests are exact. Normalization transport, the interpretation of
the conductor resultant, and the Orevkov local-orbit argument are written
proof steps. Full formal verification, an independent end-to-end proof,
external review, novelty priority, and a proof of JC(2) are not claimed.

The new scoped checks pass. Repository-wide navigation validation still
stops at pre-existing absent artifact references, also observed before the
new entries were added. The new entries and methods were checked separately,
and navigation was regenerated through the repository renderer functions.
See [`../validation.json`](../validation.json) for the exact boundary; missing
historical evidence was not reconstructed to force a green check.
