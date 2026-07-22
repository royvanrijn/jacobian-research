# Open problems for the cancellation programme

The cancellation family has separate construction, arithmetic, boundary, and
rigidity theorems.  Arithmetic classification and universality beyond the
current skeleton remain open.  The formerly open parameter-equivalence
problem is now closed by reconstruction-open faithfulness.

## 1. Closed: parameter equivalence

For fixed `(m,r)`, the common finite normalization together with its marked
affine reconstruction open recovers the selected parameter root.  The visible
target quotient remains the weighted scaling torus

\[
 (Q,R)\mapsto(uQ,u^{-m}R),
\]

and it fixes every parameter root.  It is no longer necessary to prove that
the residual target congruence kernel is trivial.  Any element of that kernel
which lifts and preserves the affine reconstruction open must preserve the two
affine factors `P=AB`, and the source identities force

\[
 A\mapsto A,\qquad y\mapsto uy,\qquad B\mapsto u^{m+1}B.
\]

It therefore fixes the reconstruction-pole residue

\[
 q=\left.\frac{B}{y^{m+1}}\right|_{A=0}.
\]

Hence the `mr` distinct geometric roots of `M_(m,r)` give `mr` distinct
stable left--right classes.  The full proof, including stabilization and the
`m=1` factor-swap exclusion, is
[the cancellation-parameter faithfulness theorem](../papers/marked-root-multiplicity/cancellation-parameter-faithfulness.tex).
Computing the target kernel itself remains an optional automorphism-group
problem, not an input to parameter faithfulness.

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
affine-reciprocal form of the same irreducibility problem is the normalized
derivative

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

A first part of the sharper fixed-row target is now proved.  The same
Newton-polygon argument gives, for density-one many `m`, a prime
`p>(log X)^10` dividing the Galois-group order and an order-`p` permutation
of type `p^t 1^f` with `f>=p-1`; see
[FIXED_R_NEWTON_RAMIFICATION.md](FIXED_R_NEWTON_RAMIFICATION.md).  It does not
yet give a Jordan cycle, because inertia can act on several root-of-unity
clusters at once.  The two precise remaining tasks are to isolate one
`p`-cycle (`t=1`) and to prove primitivity.  Their combination would imply
that the natural group contains `A_(mr)` for density-one many `m`.  Since
square discriminants have count `O_r(sqrt(X))`, this would already give the
cleaner density-one conclusion `Gal(M_(m,r))=S_(mr)`.  The alternating square
locus and the remaining reducibility exceptions can then be attacked as
separate thin sets.

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

The [endpoint-moment reduction](CONTACT_RESULTANT.md) now gives four further
uniform theorems independent of parameter irreducibility.  It proves
nonvanishing for every `m` in the columns `r=1,2,3,4`; for `r=1` it also gives
the closed resultant `((m+1)(m+2))^(-m)`.  For fixed `r`, the same method
replaces the growing degree-`mr` Sylvester problem by two degree-`r` equations
in `y=1-w` and `z=y^m`.  In the `r=3` column, Schur--Cohn separates the
degree-six endpoint eliminant from the uniform negative-binomial root disk.
In the `r=4` column, Schur--Cohn inertia isolates one exceptional conjugate
pair, and rational Rouche, angle, and Bernstein certificates prove that its
linear-subresultant `z`-branch cannot equal `y^m`.

The remaining uniform theorem, now beginning at `r>=5`, is the nonvanishing
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

The orientation-reversing reciprocal branch is now closed after its
height-one link and boundary-noncontraction marking are supplied.  Unsliced
Hensel rigidity in the whole source algebra forces the cancellation jet, a
global slice, full Stein degree one, and the two primitive rank-two end flags.
Thus the remaining cancellation-side problem is upstream: prove that an
arbitrary divisor-minimal suspension supplies that reciprocal height-one
link.  The orientation-preserving weighted chart classification remains
separate; see the
[log-geometric suspension bridge](LOG_GEOMETRY_OF_SUSPENSIONS.md).
