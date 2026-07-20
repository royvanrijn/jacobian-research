# A three-dimensional counterexample to the Jacobian conjecture

## Abstract

This repository gives an exact certificate for a polynomial map

\[
F:\mathbb A^3\longrightarrow\mathbb A^3
\]

with constant nonzero Jacobian determinant and a fiber containing three
distinct rational points. It therefore supplies a counterexample to the
classical Jacobian conjecture in dimension three. The inverse problem reduces
to a cubic polynomial, making it possible to determine the exact image, fiber
stratification, nonproperness set, normalization, and monodromy.

The same construction belongs to a weighted family controlled by the
one-variable pencil `H(W)-sW+t`. For this family the repository proves full
symmetric monodromy, a generic nodal-cuspidal discriminant theorem in every
degree, generic surjectivity from inverse degree five onward, a complete
contact-partition description of omitted values, and the associated
finite-field Chebotarev law. Every main algebraic claim has an executable exact
certificate or an explicitly identified external theorem input.

## 1. Main theorem

\[
F(x,y,z)=\left(
(1+xy)^3z+y^2(1+xy)(4+3xy),
y+3x(1+xy)^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

Exact expansion gives

\[
\det DF=-2.
\]

Nevertheless,

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)
=(-1/4,0,0).
\]

Thus `F` is everywhere étale over `C` but is not injective and
hence is not a polynomial automorphism. Appending identity coordinates gives
counterexamples in every dimension at least three. The coordinate degrees are
`(7,6,4)`.

The determinant, collision, and degrees are checked twice: once with SymPy and
once with an independent standard-library sparse-polynomial implementation.
No global geometry or family theorem is needed for this minimal certificate.

## 2. Inverse equation and reconstruction

Write a target as `(a,b,c)` and, on `x!=0`, set

\[
T=y+{1\over x}.
\]

The inverse problem becomes

\[
P(T)=cT^3-2T^2+bT-2a=0.
\]

If

\[
r=P'(T)=3cT^2-4T+b,
\]

then a simple root reconstructs the source by

\[
x={2\over r},\qquad
y=T-{r\over2},\qquad
z={5r^2\over4}-{3Tr\over2}-{cr^3\over8}.
\]

The function-field extension therefore has degree three. Its discriminant is

\[
\operatorname{Disc}_T(P)=-4Q,
\]

where

\[
Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

This model also explains the failure of properness: repeated inverse roots are
precisely the reconstruction poles, so sheets can escape to infinity without
finite ramification.

## 3. Image, fibers, and nonproperness

Define

\[
\Gamma=V(3bc-4,12a-b^2).
\]

Then

\[
F(\mathbb C^3)=\mathbb C^3\setminus\Gamma,
\qquad
S_F=V(Q),
\]

and `Gamma=Sing(V(Q))`. The complete affine fiber table is:

| Target stratum | Fiber cardinality |
|---|---:|
| `Q!=0` | 3 |
| `Q=0` away from `Gamma` | 1 |
| `Gamma` | 0 |

The proof includes the `x=0` chart, every exceptional cubic, explicit escaping
paths over all boundary strata, and a converse boundedness argument off
`V(Q)`. Root meridians generate the full monodromy group `S_3`.

See [Construction and anatomy](notes/CONSTRUCTION.md) and
[Exact image, fibers, and nonproperness](notes/IMAGE_AND_NONPROPERNESS.md).

## 4. Weighted inverse families

The three-dimensional map is part of a weighted construction determined by an
admissible primitive `H` of degree `n`. Its inverse pencil is

\[
E_{s,t}(W)=H(W)-sW+t,
\]

with `(s,t)=(BC,cAC^2)` on `C!=0`. Simple roots reconstruct finite source
points; repeated roots are reconstruction poles. The repeated-root
discriminant has normalization

\[
\nu_H(r)=\bigl(H'(r),rH'(r)-H(r)\bigr),
\]

the curve of tangent lines to the graph `Y=H(W)`.

The repository proves:

1. The generic inverse pencil is irreducible, and its geometric and arithmetic
   monodromy groups are `S_n`.
2. For generic admissible `H`, the projective discriminant is a rational
   degree-`n` curve with a smooth point at infinity, exactly `n-2` ordinary
   cusps, and

   \[
   {(n-2)(n-3)\over2}
   \]

   ordinary nodes, with no other singularities.
3. Generic weighted maps are surjective over the algebraic closure for every
   inverse degree `n>=5`.
4. The canonical family `H_d(W)=W^d(1-W)`, one-extra-root deformations, and
   repeated primitive roots have exact image and boundary theorems.

The main references are
[the weighted-seed theorem](notes/WEIGHTED_SEED_THEOREM.md),
[the generic discriminant theorem](notes/GENERIC_DISCRIMINANT_CURVE.md),
[the canonical image theorem](notes/CANONICAL_FAMILY_IMAGE.md), and
[the repeated-root boundary theorem](notes/REPEATED_ROOT_BOUNDARY.md).

## 5. Contact partitions and exceptional seeds

An omitted value on `C!=0` occurs exactly when every root of the inverse
polynomial is multiple. For a full contact partition

\[
\lambda=(\lambda_1,\ldots,\lambda_k),
\qquad \lambda_i\ge2,
\qquad \sum_i\lambda_i=n,
\]

put

\[
M_\lambda(W)=\prod_i(W-r_i)^{\lambda_i},
\qquad
\Phi_\lambda=M_\lambda(1)-M_\lambda(0)-M_\lambda'(0).
\]

The API quotients roots with equal multiplicity before elimination, saturates
all diagonals and weighted-admissibility factors, and returns the exact
coefficient-space elimination presentation. It also supports partial contact
with a residual factor.

The first exceptional parameter spaces are now explicit:

- Degree five: the `(3,2)` incidence recovers the established polynomial
  `F(R,P)` up to the scalar `5/64`.
- Degree six: the main `(2,2,2)` locus is an irreducible rational quartic
  surface; the exact `(2,2,2)` and `(3,3)` strata are disjoint and their
  closures meet exactly along the `(6)` collision boundary.
- Degree seven: the leading `(3,2,2)` locus has codimension two in normalized
  seed space. Hence nonsurjectivity is not detected by a single coefficient
  equation.
- In every degree, the exceptional locus is the union of the full-contact
  strata, with

  \[
  \dim E_\lambda=\ell(\lambda)-1,
  \qquad
  \operatorname{codim}_{\mathcal A_n}E_\lambda=n-\ell(\lambda)-2.
  \]

  The proof is a weighted-Vandermonde rank calculation on the top
  coefficients of `M_lambda`.
- Degree eight tests the proposed component poset: the maximal candidates
  `(2,2,2,2)` and `(3,3,2)` share the `(6,2)` collision boundary and have no
  exact off-diagonal intersection.
- In all degrees, merging parts defines a collision partial order with
  `E_mu` contained in the closure of `E_lambda` whenever `lambda<=mu`.
  Mason--Stothers excludes every off-collision two-omission incidence between
  distinct maximal 2/3 partitions. The maximal root hypersurfaces are
  uniformly irreducible, so the irreducible components of the exceptional-
  locus closure are indexed exactly by partitions using only twos and threes.

See [Uniform exceptional seeds](notes/UNIFORM_EXCEPTIONAL_SEEDS.md),
[contact-partition strata](notes/CONTACT_PARTITION_STRATA.md), and
[omitted-value classification](notes/OMITTED_VALUE_CLASSIFICATION.md).

## 6. Arithmetic and normal forms

At primes of good reduction, `S_n` monodromy gives the limiting fixed-point
law of a random permutation. If `D_n` is the derangement number, then the
limiting image density is

\[
1-{D_n\over n!},
\]

and targets with exactly `j` rational preimages have count

\[
{\binom{n}{j}D_{n-j}\over n!}q^3+O_H(q^{5/2}).
\]

The original cubic also has exact finite-field distributions in every
characteristic. Separately, standard reductions produce explicit
95-dimensional cubic-homogeneous and 510-dimensional Drużkowski forms with
transported rational collisions.

See [Finite-field Chebotarev](notes/FINITE_FIELD_CHEBOTAREV.md),
[exact cubic finite-field distributions](notes/FINITE_FIELD_VALUE_DISTRIBUTION.md),
[the cubic-homogeneous reduction](notes/CUBIC_HOMOGENEOUS_REDUCTION.md), and
[the cubic-linear reduction](notes/CUBIC_LINEAR_REDUCTION.md).

## 7. Reproducibility

The minimal certificate requires only the Python standard library:

```bash
python3 scripts/verify_counterexample_independent.py
```

The complete exact suite is:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
make verify
```

Useful targets are:

| Command | Scope |
|---|---|
| `make verify-minimal` | dependency-free determinant, collision, and degrees |
| `make verify-core` | minimal certificate plus cubic inverse identities |
| `make verify-geometry` | cubic image and nonproperness |
| `make verify-family` | weighted families, discriminants, contact strata, and quartic geometry |
| `make verify-normal-forms` | regenerate and verify the large normal forms |
| `make scan-weighted-seeds` | reproduce the exploratory bounded seed scan |

Julia is optional and is used only for the numerical continuation benchmark in
[NONPROPER_FIBER_BENCHMARK.md](notes/NONPROPER_FIBER_BENCHMARK.md).

## 8. Repository guide and scope

The logical status of every claim is separated from discovery provenance:

- [FACTS.md](notes/FACTS.md) records the principal exact statements.
- [IMPLEMENTATION_STATUS.md](notes/IMPLEMENTATION_STATUS.md) maps statements to
  executable support.
- [CHECKLIST.md](notes/CHECKLIST.md) distinguishes completed certificates from
  remaining independent audits.
- [PROVENANCE_AUDIT.md](notes/PROVENANCE_AUDIT.md) and
  [SOURCES.md](notes/SOURCES.md) record the available historical trail.

Additional exact analyses include the commuting inverse-Jacobian frame
([COMMUTING_FLOWS.md](notes/COMMUTING_FLOWS.md)), gradient dynamics at infinity
([GRADIENT_INFINITY.md](notes/GRADIENT_INFINITY.md)), and audited downstream
implications ([DIRECT_CONSEQUENCES.md](notes/DIRECT_CONSEQUENCES.md)).

The 627-seed scan is retained as an exploratory ledger, not as evidence for
the theorem-level claims. The foundational three-dimensional certificate does
not depend on the family theory, numerical continuation, high-dimensional
reductions, or provenance assertions.
