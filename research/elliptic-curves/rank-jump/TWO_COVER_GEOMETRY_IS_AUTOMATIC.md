# The common geometry of the 2-covers is automatic

The 48 retained production classes, including 22 exceptional-quotient
lifts, and the two small rational/Sha control classes all satisfy the
same exact trace-form identities. Their 100 signed 2-cover models are
smooth genus-one curves. Within a fixed elliptic curve, the normalized
quadric-pencil determinant does not depend on the class.

A genuine 2-cover of a smooth elliptic fibre does not gain rational
components or drop genus when it becomes rationally soluble. Those
geometric events cannot distinguish Mordell--Weil classes from Sha.
Auxiliary curves other than the 2-covers themselves can still undergo
useful degenerations; they need their own construction and proof.

The surviving issue is arithmetic: how to find rational points on several
different genus-one double covers of rational conics at once. The common
determinant and the common cubic algebra do not resolve it.

## Exact equations and determinant theorem

Let \(k\) have characteristic zero, let
\[
f(X)=X^3+aX^2+bX+c
\]
be separable, and put \(K=k[\theta]/(f)\). For
\(\beta\in K^\times\) of square norm, write
\[
\beta(z_0+z_1\theta+z_2\theta^2)^2
 =Q_{0,\beta}+\theta Q_{1,\beta}+\theta^2 Q_{2,\beta}.
\]
Our Gram convention is \(Q_{j,\beta}(z)=z^TM_j(\beta)z\), with no
factor of \(1/2\). The two signs give
\[
D_{\beta,\sigma}:\quad
Q_{2,\beta}(z)=0,\qquad Q_{1,\beta}(z)+\sigma h^2=0,
\qquad \sigma=\pm1.
\]
If \(N\beta=n^2\), taking norms gives the elliptic image
\[
x=\frac{Q_{0,\beta}}{h^2},\qquad
y=\frac{nN(z_0+z_1\theta+z_2\theta^2)}{h^3},
\]
on
\[
E_\sigma:\quad y^2=x^3+\sigma ax^2+bx+\sigma c.
\]
Thus \(\sigma=1\) is the original cubic and \(\sigma=-1\) its monic
negative twist. These are the usual labelled 2-cover equations.

Let \(H\) denote multiplication by \(\theta\), \(M_\beta\) multiplication
by \(\beta\), and \(T=H+aI\). Directly,
\[
H=\begin{pmatrix}0&0&-c\\1&0&-b\\0&1&-a\end{pmatrix},
\qquad
S=M_2(1)=
\begin{pmatrix}
0&0&1\\
0&1&-a\\
1&-a&a^2-b
\end{pmatrix}.
\]
Then
\[
\boxed{M_2(\beta)=SM_\beta,\quad
M_1(\beta)=M_2(\beta)T,\quad
\det M_2(\beta)=-N\beta.}
\]

For the four-variable Gram matrices
\[
A_\beta=\operatorname{diag}(M_2(\beta),0),\qquad
B_{\beta,\sigma}=\operatorname{diag}(M_1(\beta),\sigma),
\]
this proves
\[
\boxed{\det(\lambda A_\beta+\mu B_{\beta,\sigma})
 =-\sigma N\beta\,\mu
 \bigl[\lambda^3+2a\lambda^2\mu
 +(a^2+b)\lambda\mu^2+(ab-c)\mu^3\bigr].}
\]
After division by \(-\sigma N\beta\), no class-dependent quantity remains.

### Proof and smoothness

Coefficient extraction in a separable cubic has the trace descriptions
\[
[\,\theta^2\,]v=\operatorname{Tr}_{K/k}\frac{v}{f'(\theta)},\qquad
[\,\theta\,]v=\operatorname{Tr}_{K/k}
 \frac{(\theta+a)v}{f'(\theta)}.
\]
Multiplication is self-adjoint for these trace forms. This gives both
matrix identities. Since \(\det S=-1\) and
\(\det M_\beta=N\beta\), it gives the conic determinant. Finally
\[
\det(\lambda I+\mu T)
=N_{K/k}\bigl(\lambda+\mu(\theta+a)\bigr)
\]
is the displayed cubic factor.

The first singular member of the pencil is the rank-three cone at
\(\mu=0\); the other three have parameters
\(\lambda/\mu=-a-\theta_i\), where the \(\theta_i\) are the roots of \(f\).
All four are distinct. Their splitting field and permutation action
are determined by the original cubic, independently of \(\beta\).

The binary-quartic determinant computes the genus-one invariants of a
pair of quadrics; nonzero discriminant gives smooth genus one. See
[Fisher, *The invariants of a genus one curve*, Theorem 4.4 and §7.3](https://www.dpmms.cam.ac.uk/~taf1000/papers/g1inv.pdf).
Fisher uses Hessian matrices, twice our Gram matrices; this changes
scalars but not nonvanishing or pencil roots. Here separability of \(f\)
and \(N\beta\ne0\) prove the required nonvanishing.

Genus reduction or rational components in these models require leaving
this smooth 2-cover situation. A raw model can degenerate when its
representative ceases to be invertible, but that does not show that a
genuine 2-cover of a smooth elliptic fibre became a rational curve.
Geometrically, such a cover is always an elliptic curve.

## The first conic also fails to distinguish solubility

The conic \(Q_{2,\beta}=0\) is nonsingular because its determinant is
\(-N\beta\ne0\). If \(D_{\beta,\sigma}\) is everywhere locally soluble
over \(\mathbb Q\), projection to \(z\) supplies a point on this conic
over every completion: \(z=0\) would force \(h=0\), which is not a
projective point. Hasse--Minkowski therefore gives a rational point on
the conic.

Thus this conic is rational for every Selmer class, including nonzero
Sha classes. Parametrizing it leaves a genus-one double cover of
\(\mathbb P^1\). Finding a square there remains the solubility problem.
This generalizes the conic obstruction for the particular
[Jacobian projection family](EXPLICIT_PROJECTION_FIBRES.md) directly
to the original curves' cubic 2-cover models.

## Identical determinants do not identify the labelled covers

There is an exact criterion for a common coordinate identification.
Keeping the cubic operator \(T\) and the ordered ternary pair fixed,
there exists \(U\in\operatorname{GL}_3(k)\) such that
\[
U^TM_j(\beta)U=M_j(\gamma)\quad(j=1,2)
\]
if and only if \(\gamma/\beta\) is a square in \(K^\times\).

The two congruences imply \(U^{-1}TU=T\). The vector \(1\) is cyclic
for \(T\): the coordinates of \(1,\theta+a,(\theta+a)^2\) have determinant
one. Thus the centralizer of \(T\) is the multiplication algebra \(K\),
and \(U=M_u\) for some \(u\in K^\times\). The trace formula gives
\[
M_u^TM_j(\beta)M_u=M_j(\beta u^2).
\]
Injectivity of \(\beta\mapsto SM_\beta\) forces \(\gamma=\beta u^2\).
The converse is the same identity. Allowing a common scalar in both
congruences does not change the criterion when both norms are squares:
its cube must be a square after taking norms, so the scalar itself is
a square and can be absorbed into \(u\).

This is a **label-preserving ternary-pair criterion**, not a classification
of all projective equivalences of the four-variable models or all
isomorphisms of their underlying curves. Different rational Kummer
classes can have soluble covers whose underlying curves are all
isomorphic to the same elliptic curve. Failure of this labelled
congruence does not imply Sha.

The determinant has discarded the arithmetic squareclass, which survives
in the simultaneous trace forms. Requiring that squareclass to become
trivial would be too strong: nontrivial rational Mordell--Weil classes
also give soluble covers.

## Three paired production comparisons

The [protocol](TWO_COVER_PENCIL_GEOMETRY_PROTOCOL.json) freezes the
existing generic tested-local kernel basis and its relative lifts on
each paired control. Generic ranks and observed quotients are unchanged:

| Family / parameter | Generic rank | Observed quotient | Generic classes tested | Relative classes tested | Complete bad-place coverage |
|---|---:|---:|---:|---:|---|
| A1/MW16-05, 307/206 | 16 | 9 | 1 | 9 | yes |
| A1/MW16-05, -3158/1291 | 16 | 0 observed | 8 | 0 | no |
| A1/MW16-04, -1647/91 | 16 | 9 | 1 | 7 | no |
| A1/MW16-04, -2177/2397 | 16 | 0 observed | 8 | 0 | no |
| published R17, -2300/843 | 17 | 7 | 2 | 6 | yes |
| published R17, -1561/3133 | 17 | 0 observed | 6 | 0 | yes |

The 48 classes comprise 26 generic classes and 22 relative lifts.
Incomplete rows remain tested-place kernels, not certified full strict
spaces. These finite bases do not exhaust a full Selmer group.
Observed-zero controls retain no full curve-rank upper bound.

For every row and both signs, exact arithmetic verifies the common
operator, normalized determinant and smoothness. All original-sign
classes are rational Kummer classes by their retained point-subgroup
construction. The negative-sign covers are **not** all declared rational
or even locally soluble by this calculation. The
[complete production twist and CT results](PRODUCTION_TWIST_INCIDENCE_AND_SOLUBILITY.md)
remain the authority for those arithmetic statements.

The six normalized determinants differ with their elliptic cubics, as
they must. What is absent is any additional class-specific splitting,
rank drop or component singling out the exceptional quotient. Its
absence follows from the universal identity, not a statistical null.

### Exact rational-versus-Sha falsification

For the small control,
\[
f(X)=X^3-11X^2-14X-1,\qquad
\beta_0=\theta^2-10\theta+1,\quad
\beta_1=\theta^2-13\theta+12,
\]
both norms equal \(625\). The normalized pencil for both classes and
both signs is
\[
\mu(\lambda^3-22\lambda^2\mu+107\lambda\mu^2+155\mu^3).
\]
The original two covers are everywhere locally soluble but represent
nonzero Sha classes. The negative twist's covers are rationally soluble,
with the same first conics and explicit points
\[
(z_0,z_1,z_2,h)=(2,11/5,-1/5,1),\quad(2/5,1/5,0,1).
\]
The new verifier checks these conic points and the signed second
equation. The exact rank/CT proof is retained in
[the small rational/Sha certificate](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md).
No rank or CT result is recomputed or inferred from the pencil.

Smoothness, the rational first conic, cubic splitting field and normalized
determinant all survive a certified two-dimensional rational/Sha switch.

## Mechanism ranking and what remains

1. **Still viable — simultaneous arithmetic solubility:** a shared
   parameter cover can supply independent sections, as the split-cubic
   model proves. Its production analogue must force rational points on
   the residual genus-one covers, not merely their first conics.
2. **Established production structure — incidence:** strict cubic class
   blocks and local/CT comparisons retain information that the common
   geometric determinants discard.
3. **Excluded as distinguishing features:** genus reduction or rational
   components of a genuine 2-cover of a smooth fibre; class-specific
   pencil splitting beyond the cubic; rationality of the first conic
   for a Selmer class. The rational/Sha controls explicitly refute the
   proposed solubility implications.
4. **Missing computation:** an arithmetic relation between several
   residual double-cover equations yielding simultaneous rational
   points without inserting those points into the construction.
   Any common auxiliary curve must resolve this residual obstruction,
   rather than merely reproduce the shared determinant.
5. **For Agent 1:** no selector follows. The common determinant is
   structural bookkeeping, the residual rational-point condition is
   **solubility**, and coordinate ease is **visibility**. These automatic
   geometric properties add no **incidence** signal for new directions.

## Certificates and replay

- [Frozen class inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_two_cover_pencil_inputs_v1.json)
- [All 100 signed cover models](../../artifacts/generated-results/elliptic-curves/rank_jump_two_cover_pencil_geometry_v1.json)
- [Universal identities and independent trace verification](../../artifacts/generated-results/elliptic-curves/rank_jump_two_cover_pencil_verification_v1.json)
- [Independent class provenance and quotient accounting](../../artifacts/generated-results/elliptic-curves/rank_jump_two_cover_class_source_verification_v1.json)

The producer constructs forms by multiplication in the cubic quotient.
The verifier independently reconstructs them using traces and the inverse
of \(f'(H)\). It checks the determinant identity universally over an
eight-variable polynomial ring, cubic irreducibility at a small prime
in every case, and the explicit small conic points. Capture has a
30-second cap and an exclusive checkpoint per case. No factorization
campaign, point search or parameter sweep is involved.
The separate source verifier recomputes all products defining the 48
production classes using Sage cubic arithmetic, checks their original
independence fingerprints and certifies the 22 relative quotient lifts.

    sage -python elliptic-curves/rank-jump/two_cover_pencil_geometry.py check
    sage -python elliptic-curves/rank-jump/verify_two_cover_pencil_geometry.py check
    sage -python elliptic-curves/rank-jump/verify_two_cover_class_sources.py check

All point-derived classes are explicitly retrospective oracle data.
Only new rank-jump files, immutable certificates and this analysis index
are changed; Agent 1's work and mathematical-status entries are untouched.
