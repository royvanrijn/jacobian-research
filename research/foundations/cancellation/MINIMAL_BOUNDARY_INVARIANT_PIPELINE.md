# Operational minimal-boundary invariant pipeline

## 0. Status and scope

This note turns the eight predicates of
[the minimal-boundary classification program](MINIMAL_BOUNDARY_CLASSIFICATION.md)
into a finite certificate algorithm once the canonical normalization has
been exported as exact prime, valuation, collision, and marked-chart data.
The implementation is
[`jcsearch/minimal_boundary.py`](../jcsearch/minimal_boundary.py), with exact
family adapters in
[`jcsearch/minimal_boundary_examples.py`](../jcsearch/minimal_boundary_examples.py).

This is a **computation**, not the missing minimal-boundary theorem.  In
particular, it does not compute

\[
\operatorname{Norm}_{\mathbb A^3} k(\mathbb A^3)
\]

from an arbitrary polynomial Keller map, and it does not extract the quotient
flag or conormal marking from an unmarked normalization.  The input boundary
is therefore:

\[
\boxed{
\text{finite canonical-normalization export}
\longrightarrow
\text{eight certified decisions and a straightening attempt}.
}
\tag{0.1}
\]

The conjectural first arrow from a bare boundary-minimal map to that finite
export remains open.

## 1. Finite input

The pipeline accepts the following exact data.

1. Every height-one boundary prime has its color, ramification index,
   residue degree, different exponent, completed-incidence signature,
   quotient image, target image, and automorphism-orbit size.
2. Every declared localization link records normality, factoriality, scalar
   units, primality of its two boundary equations, the unit-lattice rank,
   the two cross-supports, its orientation exponent, and valuations of a
   finite algebra-generating set.
3. Every prime on the normalized graph has one row
   \[
   \left(
   v_V(R_{q_X}),\
   v_V(q_X^*R_\Phi),\
   v_V(F^*R_{q_Y})
   \right),
   \tag{1.1}
   \]
   together with exceptional and localization flags.  Unrecorded Rees
   primes are a separate explicit list.
4. The selected critical normalization has its parameter, genus, punctures,
   and exact unit generators.  The pipeline computes their integer valuation
   rows.
5. The conormal record gives local content vectors, conormal generator
   matrices at every recorded collision, and an exact rational residue
   coefficient.  The pipeline computes content primitivity, collision ranks,
   and the residue divisor.
6. Noncontraction is tested on exact rational residue generators in one
   normalization parameter.
7. The final optional chart record consists of marked source variables,
   target expressions, the controlled divisor and exponent, and the
   Fitting-support rows.  It contains no family label.

Family-specific calculations are confined to adapters which produce this
record.  The extractor does not import those adapters.

## 2. Algorithm

### 2.1 Intrinsic critical-prime selection

Filter the boundary primes by critical Fitting support.  A prime survives
only if its automorphism orbit has size one.  Selection succeeds precisely
when there is one such fixed prime and it is the sole critical boundary over
its quotient image.  The certificate returns its complete signature, not its
input label.

This implements `SCB` for the finite export.  It does not prove that a bare
normalization has a unique candidate.

### 2.2 Punctures and valuation rows

For \(s\) punctures, every supplied unit row must have length \(s\) and sum
zero.  Exact integer row reduction computes its rank.  The export is complete
when that rank is \(s-1\).  Together with geometric integrality, smoothness,
and genus zero, `PR<=1` passes exactly for ranks zero and one.

The output retains the rows, so `A^1` and `G_m` cannot be confused by a bare
puncture count.

### 2.3 Saturation and monotonicity

For each link, `SAT` checks all ring flags, both prime equations, a rank-one
unit lattice, singleton cross-supports, and orientation exponent
\(\epsilon=\pm1\).  The sign is read from the link data rather than chosen
from the desired chart.

For `BM`, the appropriate oriented list of algebra-generator valuations is
tested for nonnegativity.

### 2.4 Complete minimal ledger

For every graph prime \(V\), the algorithm checks

\[
v_V(R_{q_X})+v_V(q_X^*R_\Phi)-v_V(F^*R_{q_Y})=0.
\tag{2.1}
\]

It rejects:

- a nonzero signed row;
- a recorded prime which is neither exceptional, a localization boundary,
  nor present in an absolute divisor; or
- any declared unrecorded graph prime.

Thus a zero row is not by itself permission to add a spectator prime.

### 2.5 Primitive conormal and noncontraction

`PC` requires height-one saturation index one, equality of the nilradical and
conormal-generated dimensions at every collision, and a primitive residue
mark.  In rank zero the residue coefficient must have affine degree one.  In
rank one its puncture divisor must be a primitive degree-zero integer row.

`NC` differentiates the exact residue generators in the normalization
parameter.  At least one nonconstant generator gives transcendence degree
one; a constant spectator mark is rejected.

### 2.6 Coefficient straightening

The straightener does not compare against stored example formulas.

For a positive two-variable chart \((w,q)\mapsto(q,T)\), it computes the
Jacobian and tests whether

\[
\det D(q,T)=u(q-h(w)).
\tag{2.2}
\]

It then extracts \(h\), integrates it, and checks coefficientwise that

\[
T=u'\left(wq-\int h(w)\,dw\right)+g(q).
\tag{2.3}
\]

The returned \(g(q)\) is the removable target shear.

For a reciprocal chart \((P,S,Q)\), it first tests the exact controlled
Jacobian quotient.  If the second target is \(Q\), differentiation of the
third target at fixed \(P,Q\) extracts the reciprocal integral chart.  If
the second target is \(B=Q+\beta(P,S)\), changing coefficientwise from \(Q\)
to \(B\) tests

\[
C=Y(P,S)-B X(S).
\tag{2.4}
\]

This extracts \(X\) without being told its degree.

## 3. Exact regression

Run

```bash
.venv/bin/python scripts/verify_minimal_boundary_pipeline.py
```

or regenerate the checked result ledger with

```bash
.venv/bin/python scripts/verify_minimal_boundary_pipeline.py \
  --write-artifact
```

The generated ledger is
[`artifacts/generated-results/minimal_boundary_pipeline.json`](../artifacts/generated-results/minimal_boundary_pipeline.json).
The checked-in artifact has SHA-256
`6ee910a1a0244ae11c4c8d4254351897a3cb8a688e09768bc77f14b305564838`;
it uses the repository Python lock and exact characteristic-zero SymPy
arithmetic.
It records the command-independent predicate certificates for:

- weighted tangent degrees \(3,\ldots,8\);
- cancellation parameters
  \[
  (m,r)=(1,1),(2,1),(3,1),(1,2),(2,2),(1,3);
  \]
- quadratic-gauge degrees \(3,\ldots,8\);
- a chart perturbation;
- an imprimitive conormal perturbation;
- a contracted residue mark;
- a nonsaturated link; and
- a spectator critical prime and redundant zero ledger row.

The family ranges are regression representatives.  Their adapter formulas
are uniform in the parameters; the artifact is not an enumeration of an
infinite family.

As a blindness check, the verifier erases suggestive names, relabels all
primes and rows, reverses their order, and obtains the same predicate vectors
and extracted chart mechanisms.

## 4. Result

The exact outcomes are:

\[
\begin{array}{c|c|c|c}
\text{input}&\text{first seven predicates}&
\text{extracted chart}&\text{MBPkg}\\ \hline
\text{weighted}&\text{pass}&X=w&\text{pass}\\
\text{cancellation}&\text{pass}&X=S&\text{pass}\\
\text{quadratic gauge}&\text{pass}&X=S^2&\text{fail at CS}\\
\text{perturbations}&\text{localized failure}&
\text{as recorded}&\text{fail}\\
\text{spectator}&\text{SCB and LC fail}&
\text{unchanged core}&\text{fail}
\end{array}
\tag{4.1}
\]

The quadratic-gauge result is structurally important.  The coefficient
algorithm recovers the primitive quadratic mark \(S^2\), but `CS` in
Definition 1.9 of the classification program names only weighted and
cancellation diagrams.  The implementation therefore records a successful
quadratic-incidence extraction while correctly returning `CS=fail`; it does
not silently enlarge `MBPkg`.

No previously missing implication is proved by this regression.  What is
now operational is the implication from a complete finite export to eight
auditable decisions.  The smallest remaining mathematical target is still

\[
\boxed{
\text{unmarked canonical normalization}
\Longrightarrow
\text{intrinsic quotient flag and primitive marked conormal}.
}
\tag{4.2}
\]

The spectator model shows why a determinant ledger alone cannot prove this
arrow, and the quadratic-gauge row shows why chart extraction should return
a mechanism before testing membership in the two-class `CS` predicate.
