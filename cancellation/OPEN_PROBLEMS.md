# Open problems for the cancellation programme

The cancellation family has separate construction, arithmetic, boundary, and
rigidity theorems.  Unrestricted parameter equivalence, arithmetic
classification, and universality beyond the current skeleton remain open.

## 1. Finish parameter equivalence

Compute the polynomial target automorphism group preserving the labelled
boundary pair and determine which elements lift to the finite normalization.
The visible quotient is the weighted scaling torus

\[
 (Q,R)\mapsto(uQ,u^{-m}R),
\]

and it fixes every parameter root.  The decisive loose end is the residual
congruence kernel acting trivially on `P=0`.  Proving that its liftable part is
trivial would finish distinct-root left--right equivalence.

## 2. Complete the arithmetic

In standard truncated-binomial notation

\[
 P_{N,k}(x)=\sum_{j=0}^k\binom Njx^j,
\]

prove the divisibility-diagonal case

\[
 P_{N,k}\text{ irreducible whenever }
 1\le k\le N-2\quad\text{and}\quad N-k-1\mid k.
\]

Equivalently, prove irreducibility of
`P_((m+1)r+1,mr)`, hence of its reciprocal transform `M_(m,r)`, for all
`m,r>=1`.  This is a structured subfamily of the standard open conjecture for
interior truncated binomial polynomials.  Classify the natural Galois group
on the same diagonal; the discriminant has an infinite square locus, so the
expected large groups must allow both alternating and symmetric cases.  The
same polynomial is the normalized derivative

\[
 \sum_{j=0}^{mr}\binom{r+j}{r}x^j
 =\frac1{r!}\left(\frac d{dx}\right)^r
   (1+x+\cdots+x^{(m+1)r})
\]

and the Jacobi specialization `J_(mr)(x,r)`.  The derivative literature
therefore proves irreducibility for density-one many `m` on every fixed-`r`
row, but even the complete `r=1` row contains the classical open geometric-
derivative conjecture.  The complete range `mr<=30`, the full `m=1` column,
additional uniform criteria, quantitative density bounds, and known `r=1`
subfamilies are recorded in [ARITHMETIC.md](ARITHMETIC.md).

A sharper intermediate target is to combine fixed-`r` density-one
irreducibility with Newton-polygon/Jordan-cycle arguments and prove that the
natural group contains `A_(mr)` for density-one many `m`.  Square
discriminants have count `O_r(sqrt(X))` for `m<=X`, so this containment would
already imply the cleaner density-one conclusion `Gal(M_(m,r))=S_(mr)`.
The alternating square locus and the remaining reducibility exceptions can
then be attacked as separate thin sets.

Determine also the minimal field over which the full symbolic collision is
defined, rather than only the evident compositum of the parameter field and
the collision polynomial's splitting field.

## 3. Prove the all-parameter contact resultant

The full reduced/scheme/formal invariant is now active, and the generic
weighted divisorial profile has executable `(e,f)`, different, completed-DVR,
residue-factor, inertia, and degree-sum extraction.  For cancellation maps,
the generic completed contact over `P=Q=0` is also executable: it retains
`k(R)[[P,Q]]`, its critical-normalization plane `k(R)[[Y,Q]]`, the pullback
`P=(Q-Y)Y^m`, and the separate branch contributions `m^2r` and `mr`.

The full prime-intersection diagram is now executable.  With
`K=K_{m,r}` and the critical coefficient `L=L_{m,r}`, the exact certificate
`Res(K,L) != 0` proves that every geometric `K`-branch has critical contact
length `m`; distinct `K`-branches meet in the common reduced central stratum.
Since `K` is a fractional-linear transform of the parameter polynomial `M`
and `K,L` are visibly not associates, this closes the diagram in every proved
parameter-irreducibility range: all `mr<=30`, the complete `m=1` ladder, and
the uniform criteria in [ARITHMETIC.md](ARITHMETIC.md).  Direct resultants are
also checked on the regression grid `1<=m,r<=5`.

The [endpoint-moment reduction](CONTACT_RESULTANT.md) now gives three further
uniform theorems independent of parameter irreducibility.  It proves
nonvanishing for every `m` in the columns `r=1,2,3`; for `r=1` it also gives
the closed resultant `((m+1)(m+2))^(-m)`.  For fixed `r`, the same method
replaces the growing degree-`mr` Sylvester problem by two degree-`r` equations
in `y=1-w` and `z=y^m`.  In the `r=3` column, Schur--Cohn separates the
degree-six endpoint eliminant from the uniform negative-binomial root disk.

The remaining uniform theorem, now beginning at `r>=4`, is the nonvanishing
problem

\[
 \operatorname{Res}(K_{m,r},L_{m,r})\ne0.
\]

A proof would remove the last per-parameter certificate from the diagram.
Completed finite-stratum conductors, restrictions of the different, valuation
filtrations, and local monodromy are refinements of the now-known incidence
graph rather than missing prime intersections.

## 4. Audit the remaining all-parameter inputs

Boundary normalization and exhaustion have a clean-room audit.  The thick
intersection formula now has independent resultant and completed-local-ring
proofs, together with bounded symbolic and standard-library exact-arithmetic
regressions.  The most valuable independent regressions still missing are:

1. the universal cancellation and reconstruction proof;
2. the all-parameter monodromy computation.

Positive invariance regressions should also conjugate a known map by
complicated source and target automorphisms, change primitive elements and
compactification coordinates, and stabilize it.  The reconstructed boundary
object should remain canonically isomorphic.

## 5. Go beyond the skeleton

[RIGIDITY.md](RIGIDITY.md) exhausts additional monomial weights, finitely many
normalized resolvent factors, and arbitrary target-dependent polynomial
derivatives while retaining one triangular reconstruction variable.  A
genuinely broader classification must change that skeleton—for example by
adding another source function, source variable, or independent inverse
variable—and then rebuild the cancellation and boundary analyses.

## 6. Classify controlled-boundary suspensions

The weighted and cancellation constructions now have explicit but distinct
plane cores:

\[
 (W,s)\longmapsto(s,Ws-H(W)),
 \qquad J=H'(W)-s,
\]

and, in a family over `P`,

\[
 (s,Q)\longmapsto
 \left(Q,C\int_0^s\{1-t(Q-Pt)^m\}^r\,dt\right),
 \qquad J=-C D^r.
\]

Classify plane maps whose Jacobian is a controlled boundary power and
determine which admit a polynomial or birational suspension to a
constant-Jacobian map in one higher dimension.  The first concrete questions
are whether every one-boundary suspension is equivalent to one of these two
normal forms, and which divisor, monodromy, and reconstruction data obstruct
such an equivalence.  The precise comparison theorem and determinant ledgers
are in the
[tangent-map core note](../verified/TANGENT_MAP_CORE.md).  A two-sided
suspension-square formalism, the simple-section plane-core normal form, and
the first independent two-boundary triangular obstruction are developed in
[the controlled-boundary suspension note](CONTROLLED_BOUNDARY_SUSPENSIONS.md).
