The map was announced by Levent Alpöge and Fable; this repository
independently verifies and studies it.

# A three-dimensional counterexample to the Jacobian conjecture

## Abstract

This repository gives an exact certificate for a polynomial map

\[
F:\mathbb A^3\longrightarrow\mathbb A^3
\]

with constant nonzero Jacobian determinant and a fiber containing three
distinct rational points. It therefore supplies a counterexample to the
classical Jacobian conjecture in dimension three. ([Claim C01](CLAIMS.md#c01))
The inverse problem reduces to a cubic polynomial. ([Claim C02](CLAIMS.md#c02))
Its exact image, fiber stratification, nonproperness set, and monodromy are
also determined. ([Claim C03](CLAIMS.md#c03))

The same construction belongs to a weighted family controlled by the
one-variable pencil `H(W)-sW+t`. For this family the repository proves full
symmetric monodromy, a generic nodal-cuspidal discriminant theorem in every
degree, generic surjectivity from inverse degree five onward, a complete
contact-partition and irreducible-component theorem for exceptional seeds, and
the associated finite-field Chebotarev law.

### Proof status

The foundational three-dimensional counterexample has direct executable
certificates and a separately authored Lean 4 formalization by Dean Cureton;
see [External Lean formalization of C01](notes/LEAN_C01.md).  The all-degree
weighted-family results are conventional
mathematical theorems proved by uniform written arguments: irreducibility of
`H(W)-sW+t`, birational normalization of the discriminant, transitivity plus
local transpositions for `S_n`, Mason--Stothers separation, weighted
Vandermonde dimension calculations, and tangent-chord collision order.

The accompanying scripts verify the algebraic identities used in those
arguments and test exact examples in bounded degrees (typically through
degree eight or ten).  Those computations are regression tests and proof
audits; they are not formal machine proofs of the quantified all-degree
statements.  External theorem inputs are identified where they occur.

## 1. Main theorem

\[
F(x,y,z)=\left(
(1+xy)^3z+y^2(1+xy)(4+3xy),
y+3x(1+xy)^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

The source is the space of binary cubics with a marked simple projective
root; the Keller map forgets the marked root.  More precisely, `(x,y,z)`
marks `[1+xy:x]` on

\[
cU^3-2U^2V+bUV^2-2aV^3=0,
\]

and this identifies the source with the simple-root locus of the finite
degree-three incidence cover.  This is the main conceptual explanation of the
determinant, generic degree three, root at infinity, `3/1/0` fiber table, and
omitted triple-root locus; see
[The cubic map as a marked-root space](notes/MARKED_ROOT_MODEL.md).

Here is a short structural proof of the determinant.  On `x!=0`, use the
inverse-chart coordinates

\[
t=y+{1\over x},\qquad r={2\over x},\qquad c=F_3(x,y,z).
\]

For a target `(a,b,c)`, the inverse cubic and its derivative give

\[
a=t^2+{rt\over2}-ct^3,qquad
b=r+4t-3ct^2.
\]

Therefore

\[
\det {\partial(a,b,c)\over\partial(t,r,c)}
=\det\begin{pmatrix}
2t+r/2-3ct^2&t/2&-t^3\\
4-6ct&1&-3t^2\\
0&0&1
\end{pmatrix}
={r\over2}.
\]

Directly from the source chart,

\[
\det {\partial(t,r,c)\over\partial(x,y,z)}
=\det\begin{pmatrix}
-1/x^2&1&0\\
-2/x^2&0&0\\
2-6xy-3x^2z&-3x^2&-x^3
\end{pmatrix}
=-2x.
\]

Since `r=2/x`, the chain rule gives

\[
\det DF={r\over2}(-2x)=-2
\]

on `x!=0`, and hence everywhere by polynomial identity.

Nevertheless,

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)
=(-1/4,0,0).
\]

Thus `F` is everywhere étale over `C` but is not injective and
hence is not a polynomial automorphism. Appending identity coordinates gives
counterexamples in every dimension at least three. The coordinate degrees are
`(7,6,4)`. ([Claim C01](CLAIMS.md#c01))

The determinant also has two executable checks: one with SymPy and one with an
independent standard-library sparse-polynomial implementation.  No global
geometry or family theorem is needed for this minimal certificate.

Dean Cureton has also formalized this determinant, the explicit collision,
and a determinant-one rescaling in Lean 4.  The pinned external replication
target `make verify-lean-c01` builds his proof without vendoring its unlicensed
source; its exact theorem scope and attribution are recorded in
[the Lean C01 audit](notes/LEAN_C01.md).

## 2. Marked-root inverse model and reconstruction

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

([Claim C02](CLAIMS.md#c02))

Projectively, `P` is the binary cubic above.  Its root `[1:0]` when `c=0`
reconstructs regularly to `(0,b,a-4b^2)`, so it is exactly the `x=0` chart,
not an escaping sheet.  The marked-root theorem identifies all of `A^3` with
the simple-root locus, including this infinity chart.  Repeated marked roots
form the ramification divisor of the ambient finite incidence cover; that
divisor is removed from the source.  It is not a branch divisor of `F`, which
is everywhere etale.

This model also explains the failure of properness: repeated finite inverse
roots are reconstruction poles, so sheets can escape to infinity without
finite ramification of the affine Keller map.

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

These counts are simply the numbers of simple projective roots of the binary
cubic: three distinct roots, a double root plus a simple root, or a triple
root.  In the last case coefficient comparison gives
`3bc=4, 27ac^2=4`, equivalently `3bc=4, 12a=b^2`.

The proof includes the `x=0` chart, every exceptional cubic, explicit escaping
paths over all boundary strata, and a converse boundedness argument off
`V(Q)`. Root meridians generate the full monodromy group `S_3`.
([Claim C03](CLAIMS.md#c03))

See [Construction and anatomy](notes/CONSTRUCTION.md) and
[Exact image, fibers, and nonproperness](notes/IMAGE_AND_NONPROPERNESS.md).
The latter remains the proof of the exact image and nonproperness statements;
the [marked-root model](notes/MARKED_ROOT_MODEL.md) packages their fiber
geometry.

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
([Claim C04](CLAIMS.md#c04))

Globally, the weighted source is the regular-reconstruction open in the
normalization of the finite marked-root incidence

\[
H(W)-BCW+cAC^2=0.
\]

Over `C!=0` this is exactly the simple-root locus. Over `C=0`, normalization
is essential: finite `gamma=0` source points can lie over the multiple root
`W=0`, while an additional simple root of `H` can still be a reconstruction
pole. See [The weighted family as a normalized marked-root
space](notes/WEIGHTED_MARKED_ROOT_MODEL.md). ([Claim C04](CLAIMS.md#c04))

The repository proves:

1. The generic inverse pencil is irreducible, and its geometric and arithmetic
   monodromy groups are `S_n`. ([Claim C04](CLAIMS.md#c04))
2. For generic admissible `H`, the projective discriminant is a rational
   degree-`n` curve with a smooth point at infinity, exactly `n-2` ordinary
   cusps, and

   \[
   {(n-2)(n-3)\over2}
   \]

   ordinary nodes, with no other singularities. ([Claim C05](CLAIMS.md#c05))
3. Generic weighted maps are surjective over the algebraic closure for every
   inverse degree `n>=5`. ([Claim C06](CLAIMS.md#c06))
4. The canonical family `H_d(W)=W^d(1-W)`, one-extra-root deformations, and
   repeated primitive roots have exact image and boundary theorems.
   ([Claim C06](CLAIMS.md#c06))

The normalized graph compactification has one universal discriminant
dicritical divisor and explicit Kummer divisors over `C=0` indexed by the
primitive-root multiplicities of `H`. Their images give the nonproperness set,
and an omitted value is exactly a fiber supported entirely on the
discriminant divisor. These are precisely the polar divisors removed from the
normalized marked-root incidence to recover the affine source.
([Claim C16](CLAIMS.md#c16))

The main references are
[the weighted-seed theorem](notes/WEIGHTED_SEED_THEOREM.md),
[the generic discriminant theorem](notes/GENERIC_DISCRIMINANT_CURVE.md),
[the canonical image theorem](notes/CANONICAL_FAMILY_IMAGE.md), and
[the repeated-root boundary theorem](notes/REPEATED_ROOT_BOUNDARY.md), and
[the dicritical compactification theorem](notes/DICRITICAL_COMPACTIFICATION.md).
The incidence-space formulation is in
[the weighted marked-root theorem](notes/WEIGHTED_MARKED_ROOT_MODEL.md).

## 5. Contact partitions and exceptional seeds

An omitted value on `C!=0` occurs exactly when every root of the inverse
polynomial is multiple. Let `A_n` be the normalized admissible seed space;
it has dimension `n-3`. For a full contact partition

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

Let `N_n` be the nonsurjective seed locus and `E_lambda` the image of the
exact full-contact incidence of type `lambda`. The complete incidence
classification is

\[
\mathcal N_n=
\bigsqcup_{\substack{\lambda\vdash n\\\lambda_i\ge2}}\mathcal E_\lambda,
\qquad
\dim\mathcal E_\lambda=\ell(\lambda)-1,
\qquad
\operatorname{codim}_{\mathcal A_n}\mathcal E_\lambda
=n-\ell(\lambda)-2.
\]

([Claim C08](CLAIMS.md#c08))

### Unique omitted-value theorem

Every normalized admissible seed has at most one omitted inverse-pencil value
`(s,t)`. Indeed, two distinct omitted values would give monic full-contact
polynomials `P,Q` with `P-Q` nonzero affine-linear. Mason--Stothers excludes
all support-length cases except

\[
n\text{ even},
\qquad
P=A^2,
\qquad Q=B^2,
\qquad \deg A=\deg B=n/2.
\]

But then

\[
P-Q=(A-B)(A+B)
\]

has degree at least `n/2>=2`, again a contradiction. Hence the displayed
union is genuinely disjoint: the exact partition of the unique omitted
polynomial assigns every nonsurjective seed to one and only one stratum.
([Claim C07](CLAIMS.md#c07))

### Main theorem

For every pair `(a,b)` of nonnegative integers satisfying `2a+3b=n`, let

\[
\mathcal C_{a,b}=\overline{\mathcal E_{2^a3^b}}.
\]

Equivalently, `C_(a,b)` is the closure of the stratum whose full-contact
polynomial has the form

\[
M=Q^2R^3,
\qquad \deg Q=a,
\qquad \deg R=b.
\]

Write `C_lambda=C_(a,b)` for `lambda=2^a3^b`. The result has three parts.

**Component classification.**

\[
\boxed{
\operatorname{Irr}\bigl(\overline{\mathcal N_n}\bigr)
=\left\{\mathcal C_\lambda:
\lambda\vdash n,\ \lambda_i\in\{2,3\}\right\}.}
\]

**Closure order.** Define `lambda<=nu` when `nu` is obtained from `lambda` by
merging parts. For every full-contact partition `nu`,

\[
\boxed{
\mathcal E_\nu\subseteq\mathcal C_\lambda
\quad\Longleftrightarrow\quad
\lambda\preceq\nu.}
\]

**Intersection formula.** For maximal 2/3 partitions `lambda` and `mu`,
set-theoretically inside `A_n`,

\[
\boxed{
\mathcal C_\lambda\cap\mathcal C_\mu
=\bigcup_{\substack{\lambda\preceq\nu\\\mu\preceq\nu}}
\mathcal E_\nu.}
\]

In short: the irreducible components are indexed by 2/3 partitions, their
boundary strata by arbitrary partitions with parts at least two, and every
component intersection is determined by common coarsening. The theorem is
set-theoretic; scheme-theoretic intersections, embedded components, and
intersection multiplicities are separate questions.
([Claim C10](CLAIMS.md#c10))

### Dimension theorem

The root support has dimension `a+b`, and endpoint normalization imposes one
equation. Therefore

\[
\dim\mathcal C_{a,b}=a+b-1.
\]

Since `dim A_n=n-3`,

\[
\operatorname{codim}_{\mathcal A_n}\mathcal C_{a,b}
=n-a-b-2=a+2b-2.
\]

Maximizing `a+b` means minimizing the number of triple parts. Hence the full
nonsurjective locus has

\[
\boxed{
\operatorname{codim}_{\mathcal A_n}\overline{\mathcal N_n}
=\left\lceil{n\over2}\right\rceil-2.}
\]

There is a unique top-dimensional component: it is indexed by all double
parts when `n` is even, and by one triple part with all remaining parts double
when `n` is odd. Whenever further representations `2a+3b=n` exist, they give
genuine lower-dimensional components; the exceptional locus is then not
equidimensional.
([Claim C08](CLAIMS.md#c08))

The codimension formula also explains the sharp degree transition in the
image theorem. The exceptional closure has codimension zero in degrees three
and four, codimension one in degree five, and increasing codimension
thereafter. Thus generic surjectivity begins exactly at inverse degree five
for structural, rather than case-specific, reasons.

### Why twos and threes?

The answer is local and additive. Omission forbids simple roots, so every
allowed contact order lies in the numerical semigroup

\[
S_2=\{2,3,4,\ldots\}=\langle2,3\rangle.
\]

Collision adds contact orders. The only indecomposable nonzero elements of
`S_2` are two and three: neither can split without introducing a forbidden
simple root, while every `m>=4` splits as `2+(m-2)`. Thus the 2/3 partitions
are the atomic collision types. The dimension formula turns these atoms into
maximal component closures; they do not arise from a low-degree pattern.

Mason--Stothers has a different role. For a full-contact partition put

\[
\epsilon(\lambda)=\sum_i(\lambda_i-2)=n-2\ell(\lambda).
\]

The affine-difference incidence for two degree-`n` types would require

\[
n\le\ell(\lambda)+\ell(\mu)
=n-{\epsilon(\lambda)+\epsilon(\mu)\over2}.
\]

For distinct partitions the total excess is positive, since the all-double
partition is the unique zero-excess type. Hence Mason excludes off-collision
coexistence for every pair of distinct full-contact partitions, not only for
maximal 2/3 types.

More generally, if reconstruction required multiplicity at least `r`, the
primitive contact orders would be `r,...,2r-1`. The specific pair `{2,3}`
comes from the present threshold `r=2`; the atom principle is the universal
part.
([Claim C09](CLAIMS.md#c09))

### Universal equation and irreducibility

For a maximal partition `lambda=2^a3^b`, write `M=Q^2R^3` as above and set

\[
x=Q(0),\quad u=Q'(0),\quad X=Q(1),
\qquad
y=R(0),\quad v=R'(0),\quad Y=R(1),
\]

then the normalized incidence equation is universally

\[
\boxed{
\Phi_{2^a3^b}
=X^2Y^3-x^2y^3-2xuy^3-3x^2y^2v.}
\]

Its irreducibility proof is uniform:

- If `b>=3`, the endpoint coordinates of `R` are independent and `Phi` is
  primitive linear in `v`.
- If `a>=3`, `Phi` is primitive quadratic in the independent coordinate `X`;
  the odd `x`-valuation of its constant term makes the discriminant
  nonsquare.
- Only seven endpoint-rank cases lie outside these stable regimes, and each
  has an explicit linear or nonsquare-discriminant certificate.

([Claim C10](CLAIMS.md#c10))

### Normalization theorem

Retain the collision diagonals in `Phi_(2^a3^b)=0`, but impose `D!=0` and the
weighted admissibility condition.  The resulting quotient-coordinate
hypersurface `tilde(C)_(a,b)` is smooth.  In the stable ranges this follows
from

\[
{\partial\Phi\over\partial u}=-2xy^3
\qquad\text{or}\qquad
{\partial\Phi\over\partial v}=-3x^2y^2,
\]

because `x=Q(0)` and `y=R(0)` cannot vanish when `Phi=0` and `D!=0`.  The
same seven endpoint-rank cases saturate to the empty singular locus exactly.

Top-coefficient recovery makes the natural seed map finite, while an exact
2/3 polynomial uniquely determines `Q` and `R`, making it generically
one-to-one.  Therefore

\[
\boxed{
\widetilde{\mathcal C}_{a,b}\longrightarrow\mathcal C_{a,b}
\text{ is the normalization morphism}.}
\]

Over the generic point of the exact collision stratum of type
`(m_1,...,m_k)`, the geometric branches of the component are counted by

\[
[U^aV^b]\prod_{\rho=1}^k
\left(\sum_{2i+3j=m_\rho}U^iV^j\right).
\]

The entries label distinct generic roots, so equal multiplicities are not
quotiented by permutation.  Thus every component has an explicit smooth
normalization and its generic collision branches are combinatorial.  The
first geometric two-branch example occurs in degree twelve for
`(a,b)=(3,2)` over `E_(6,6)`.  Degree fourteen is only the first ambiguity
after additionally quotienting abstract allocation data by permutations of
equal-multiplicity roots.
([Claim C11](CLAIMS.md#c11))

This is deliberately a generic-stratum statement: at more-special points,
branches may specialize or ramify further.  Completed-local branch structure,
scheme-theoretic ramification, and conductor multiplicities remain separate
questions in general.  The first two-branch case is now explicit: generically
along `E_(6,6)`, the completed local ring of `C_(3,2)` is the fiber product of
two regular fourfold branches over
`k[[t,epsilon,eta]]/(epsilon^2,eta^2)`.  Its conductor has transverse
colength four, so the branches are quadratically tangent rather than nodal.
See [the generic degree-twelve local singularity](notes/DEGREE12_LOCAL_SINGULARITY.md).
([Claim C12](CLAIMS.md#c12))

The first higher-transfer block is also explicit. For two simultaneous
`3Q <-> 2R` transfers at one root, even allowing the normalized polynomials to
differ by an affine term, the factorization scheme is finite flat of rank four
over the reduced quadratic-factor space, with coincident-root fiber
`k[X,Y]/(X^3,XY,Y^2)`. ([Claim C17](CLAIMS.md#c17))

The three-transfer block is finite flat of rank eight over the reduced monic
cubic-factor space.  Its coincident-root fiber has Hilbert function
`(1,3,3,1)` and socle dimension two, and affine difference again adds no
points. ([Claim C19](CLAIMS.md#c19))

The four-transfer block is finite flat of rank sixteen over the reduced monic
quartic-factor space, with affine difference equal to strong equality.  Its
coincident-root fiber has Hilbert function `(1,4,6,4,1)` and socle dimension
four.  This verifies the predicted rank `2^k` through `k=4`, while showing
that the collision algebras are already far from Gorenstein.
([Claim C20](CLAIMS.md#c20))

For arbitrary pairs of allocations, strong equality decomposes formally as a
Hensel tensor product of the rootwise transfer blocks and a smooth common
factor.  The actual normalized-branch intersection is the same product once a
single global affine-difference length bound is verified.  In the first global
higher-transfer case `(2,-2)`, that bound is sixteen and the completed
intersection is `Z_2 completed-tensor Z_2`, with Hilbert function
`(1,4,6,4,1)`.  The first mixed case `(2,-1,-1)` also has length sixteen and
equals `Z_2 completed-tensor Z_1 completed-tensor Z_1`.
([Claim C18](CLAIMS.md#c18))

### Proof architecture for closures and intersections

Collision gives the forward closure implication:

\[
\lambda\preceq\mu
\quad\Longrightarrow\quad
\mathcal E_\mu\subseteq\overline{\mathcal E_\lambda}.
\]

Every partition with parts at least two is therefore contained in the closure
of a maximal 2/3 stratum obtained by splitting its larger parts. Conversely,
for two distinct maximal components, every intersection point inside `A_n`
comes from a root collision. There are no off-collision intersections
supporting two distinct omitted values of different maximal types.

Three facts make this conclusion global: tangent-chord deformation realizes
every permitted collision closure; weighted Newton sums make recovery of the
omitted roots from the top coefficients finite, excluding roots escaping at
the boundary; and polynomial Mason--Stothers excludes the remaining
off-collision two-omission incidence for every pair of distinct contact
types. These facts prove both the reverse closure implication and the
common-coarsening intersection formula.
([Claim C10](CLAIMS.md#c10))

### Number of components

The number of components is the number of representations `n=2a+3b`:

\[
c_n=\#\left\{b:0\le b\le\left\lfloor n/3\right\rfloor,
\ b\equiv n\pmod2\right\}.
\]

If `m=floor(n/3)`, this is

\[
c_n=
\begin{cases}
\left\lfloor m/2\right\rfloor+1,&n\text{ even},\\
\left\lfloor(m+1)/2\right\rfloor,&n\text{ odd}.
\end{cases}
\]

Equivalently,

\[
\sum_{n\ge0}c_nz^n={1\over(1-z^2)(1-z^3)}.
\]

Together, these statements give completeness, irreducibility, and
distinctness of every component in all degrees.
([Claim C10](CLAIMS.md#c10))

The low-degree calculations remain useful regressions: degree five recovers
the established `(3,2)` polynomial `F(R,P)`; degree six has the two components
`(2,2,2)` and `(3,3)` meeting along `(6)`; degree seven first exhibits genuine
codimension two; and degree eight has components `(2,2,2,2)` and `(3,3,2)`
meeting along the predicted collision boundary `(6,2)`.

See [Uniform exceptional seeds](notes/UNIFORM_EXCEPTIONAL_SEEDS.md),
[the contact-atom principle](notes/CONTACT_ATOM_PRINCIPLE.md),
[the unique omitted-value theorem](notes/UNIQUE_OMITTED_VALUE.md),
[component normalization](notes/COMPONENT_NORMALIZATION.md),
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

([Claim C13](CLAIMS.md#c13))

The original cubic also has exact finite-field distributions in every
characteristic. Separately, standard reductions produce explicit
95-dimensional cubic-homogeneous and 510-dimensional Drużkowski forms with
transported rational collisions.
([Claim C15](CLAIMS.md#c15))

The explicit quartic weighted model has determinant `-6`, two cusps and one
node on its discriminant, full `S_4` monodromy, and a completely classified
image and nonproperness set. ([Claim C14](CLAIMS.md#c14))

See [Finite-field Chebotarev](notes/FINITE_FIELD_CHEBOTAREV.md),
[exact cubic finite-field distributions](notes/FINITE_FIELD_VALUE_DISTRIBUTION.md),
[the cubic-homogeneous reduction](notes/CUBIC_HOMOGENEOUS_REDUCTION.md), and
[the cubic-linear reduction](notes/CUBIC_LINEAR_REDUCTION.md).

## 7. Reproducibility

The minimal certificate requires only the Python standard library:

```bash
python3 scripts/verify_counterexample_independent.py
```

The routine executable identity and regression suite is:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
make verify
```

Python is pinned in [`.python-version`](.python-version), and every Python
dependency is pinned in [`requirements.txt`](requirements.txt).  To retain a
local environment record and full verification log, run:

```bash
make verify-logged
```

This writes `artifacts/verification/environment.txt` and
`artifacts/verification/verify.log`.  CI runs the four evidence classes in
separate jobs and archives the corresponding logs for every commit.

Useful targets are:

| Command | Scope |
|---|---|
| `make verify-minimal` | dependency-free determinant, collision, and degrees |
| `make verify-core` | foundational map, marked-root model, cubic inverse, exact image, fibers, and nonproperness |
| `make verify-theorems` | exact identities and local models used by the written uniform proofs |
| `make verify-regressions` | bounded-degree family, discriminant, Chebotarev, and quartic regressions |
| `make verify-derived` | regenerate and verify the large stable-equivalent normal forms |
| `make verify-family` | backward-compatible alias for theorem identities plus regressions |
| `make verify-normal-forms` | backward-compatible normal-form target |
| `make scan-weighted-seeds` | reproduce the exploratory bounded seed scan |

Julia is optional and is used only for the numerical continuation benchmark in
[NONPROPER_FIBER_BENCHMARK.md](notes/NONPROPER_FIBER_BENCHMARK.md).

## 8. Repository guide and scope

The logical status of every claim is separated from discovery provenance:

- [CLAIMS.md](CLAIMS.md) assigns every major README theorem a unique claim ID
  and records its written proof, executable support, regression role, external
  dependencies, and independent-audit status.
- [FACTS.md](notes/FACTS.md) records the principal exact statements.
- [IMPLEMENTATION_STATUS.md](notes/IMPLEMENTATION_STATUS.md) distinguishes
  uniform written proofs from executable identities and finite-degree
  regressions.
- [CHECKLIST.md](notes/CHECKLIST.md) distinguishes completed certificates from
  remaining independent audits.
- [PROVENANCE_AUDIT.md](notes/PROVENANCE_AUDIT.md) and
  [SOURCES.md](notes/SOURCES.md) record the available historical trail.

Additional exact analyses include the commuting inverse-Jacobian frame
([COMMUTING_FLOWS.md](notes/COMMUTING_FLOWS.md)), gradient dynamics at infinity
([GRADIENT_INFINITY.md](notes/GRADIENT_INFINITY.md)), and audited downstream
implications ([DIRECT_CONSEQUENCES.md](notes/DIRECT_CONSEQUENCES.md)).

The 627-seed scan is retained as an exploratory ledger, not as evidence for
the theorem-level claims.  Likewise, passing `make verify` does not by itself
prove the all-degree family theorems; those depend on the cited written
uniform arguments.  The foundational three-dimensional certificate does not
depend on the family theory, numerical continuation, high-dimensional
reductions, or provenance assertions.
