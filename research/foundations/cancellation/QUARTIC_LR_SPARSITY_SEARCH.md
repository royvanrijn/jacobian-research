# Bounded left--right sparsity search for the ungraded quartic

## Status

This note records an exact **bounded computation**, not an absolute sparsity
theorem.  The starting map is the rational-root determinant-one quartic

\[
 t=1+xy,\qquad q=t^2z-\frac47y^2(1+3t),
\]

\[
 \left(
 -\frac12tq,\;
 y-\frac{21}{4}xq+3t^2x^2q^4,\;
 x(5-3t)+\frac74x^3z-\frac32(xq)^4
 \right),                                               \tag{1}
\]

with expanded support \((7,51,38)\) and degrees \((7,26,24)\).
The search did **not** find a smaller representative.  It did produce:

1. exact negative evidence for three natural classes of nontrivial
   polynomial left--right transformations;
2. a much smaller standalone affine-linear symmetry certificate for (1);
3. a rational left--right scaling which improves expanded coefficient and
   collision height while retaining the same support;
4. a two-move circuit search whose best nontrivial structured candidate is
   the one-term near miss \((7,51,39)\).

The affine-linear theorem is stated and proved in
[the canonical linear-torus note](NO_LINEAR_TORUS_COUNTEREXAMPLE.md#6-standalone-sparse-affine-linear-record).
The nonlinear intrinsic theorem remains separate in
[the algebraic-torus note](NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md).

## 1. Exact search space

All calculations are over \(\mathbb Q\), with expansions performed exactly.
The bounded search contains four parts.

### Elementary source shears

For each source coordinate \(v_i\), every monomial \(m\) of total degree at
most two in the other two source coordinates was tested in

\[
 v_i\longmapsto v_i+c\,m,\qquad c\in\mathbb Q.           \tag{2}
\]

There are 15 monomial directions.  For each direction, the checker constructs
the coefficient polynomials in \(c\), enumerates every rational root of every
coefficient, and tests all resulting exceptional support patterns.  Thus the
one-parameter search is exact; it is not a bounded-height sampling of \(c\).

### Elementary target shears

For each target coordinate, every monomial of total degree at most two in
the other two target coordinates was tested in the analogous triangular
target automorphism.  Again there are 15 directions.  Since the transformed
coordinate is affine-linear in \(c\), all possible support drops occur among
the exact coefficient ratios enumerated by the checker.

### Structured source-jet shears

Write \(u=xy\).  A source shear

\[
 z\longmapsto z+y^2p(u)
\]

replaces the essential boundary part

\[
 R_0(u)=-\frac47(4+3u)
\]

of \(q\) by \(R(u)=R_0(u)+(1+u)^2p(u)\).  The search enumerates every
two-term polynomial

\[
 R(u)=Au^m+Bu^n,\qquad 0\le m<n\le12,                  \tag{3}
\]

with the same value and first derivative as \(R_0\) at \(u=-1\).  These two
Hermite conditions are exactly the condition that \(p\) be polynomial.
There are 76 candidates.

### Rational diagonal left--right scalings

The two stable-moduli scalings were enumerated with positive reduced
rationals \(\alpha,\beta\) whose numerators and denominators are at most 16.
This gives 25281 exact pairs.  Two heights were recorded:

* the largest naive reduced height of an expanded coefficient;
* the largest naive reduced height of a coordinate in the four-point
  collision.

This finite grid is an optimization experiment, not an optimality theorem.

### Two-move circuits and three-term jets

A second computation allows intermediate support growth.  It tests 330
circuits of the form

\[
 v_i\longmapsto v_i+c\,m,\qquad
 F_j\longmapsto F_j+dF_k,                               \tag{4}
\]

where \(m\) has degree at most two, the nonzero reduced rational \(c\) has
numerator and denominator at most four, and \(d\) runs through every exact
coefficient ratio at which the target cleanup can cancel a term.

Independently, it exhausts all 286 three-term boundary jets

\[
 R(u)=Au^a+Bu^b+Cu^c,\qquad 0\le a<b<c\le12.            \tag{5}
\]

The two Hermite equations at \(u=-1\) solve \(A,B\) linearly in \(C\).
Every rational root of every resulting coefficient polynomial in \(C\) is
then tested.  This is exact in the free parameter and includes
multi-monomial \(z\)-shears which decompose into commuting elementary source
moves.

## 2. Outcome

For every elementary source direction (2), the unique best exceptional
parameter was \(c=0\); no nonzero parameter preserved or reduced total
support.  The same holds for all 15 elementary target directions.

Among the 76 structured jet shears (3), the unique minimum is the original
pair \((m,n)=(0,1)\), with support \((7,51,38)\).  The closest nontrivial
candidate has support \((7,51,40)\) and degrees \((10,38,36)\).

The balanced scaling optimum on the declared rational grid is

\[
 \alpha=\frac14,\qquad\beta=\frac{12}{5}.
\]

It retains support \((7,51,38)\), lowers maximum expanded coefficient height
from \(2248704\) to \(21875\), and lowers collision-coordinate height from
\(24820\) to \(19856\).  The coefficient-only optimum is

\[
 \alpha=\frac7{12},\qquad\beta=6,
\]

with coefficient height \(4648\) but collision height \(49640\).  The
balanced representative and its collision are displayed in
[the intrinsic note](NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md#a-sparser-rational-representative).

No circuit in (4) improves total support 96.  The best source-to-target
candidate has intermediate and final support \((8,61,49)\), total 118; its
optimal target cleanup is the identity.  Some larger intermediates do admit
nonzero cleanup, but none approaches the starting record.

The exact three-term search (5) finds a unique smallest nonidentity boundary
polynomial, represented redundantly by several exponent triples, with support

\[
\boxed{(7,51,39)}
\]

at

\[
 R(u)=R_0(u)+\frac{16}{7}(1+u)^2.
\]

Equivalently, this is the elementary source shear

\[
 z\longmapsto z+\frac{16}{7}y^2.                        \tag{6}
\]

Thus the enlarged structured family comes within one term but does not
improve (1).  From this near miss, every monomial target shear through degree
three fails to improve support; its best nonidentity target move has total
108.  Every elementary second source shear through degree two is also
tested at all rational exceptional parameters.  The literal inverse of (6)
returns to 96, while the best noninverse move has total 113.
Finally, adjoining a second structured monomial \(bu^n\) to (6), for every
\(1\le n\le12\), is checked at every rational exceptional \(b\).  No such
move improves 97; the best nonzero continuation has support \((7,51,40)\),
total 98.

The separate affine-linear calculation for the original sparse
representative extracts 24 primitive coefficient rows with only 46 nonzero
entries and determinant \(10\).  This is the checkable weak record requested
before the full boundary theorem: determinant one, a four-point collision,
and one small square determinant.

## 3. What the search does not exclude

The computation does not cover:

* compositions of elementary shears whose intermediate support increases;
* general two-source-move circuits away from the near miss (6);
* source-to-target circuits outside the rational parameter box in (4);
* elementary shears of source degree at least three outside the structured
  \(t\)-adic family;
* target shears of degree at least three, except for the near-miss exit audit;
* wild polynomial automorphisms;
* alternate polynomial algebraizations of the marked-line incidence;
* another genuinely ungraded left--right class.

In particular, the unchanged minimum 96 is bounded search evidence only.
The proved minimum remains the family-relative statement \(g_2=0\) inside
the displayed quadratic-gauge normal form.

The most focused next search is a two-move circuit search which permits a
temporary support increase, followed by a bounded alternate-lift search for
cross-cancellations among several coefficient weights.  The one-move and
single-jet routes recorded here are exhausted only in their stated bounds.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/search_quartic_lr_sparsity.py \
  --source-degree 2 \
  --target-degree 2 \
  --tadic-max-exponent 12 \
  --scaling-bound 16 \
  --output artifacts/generated-results/quartic_lr_sparsity_search.json
```

The generated record has SHA-256

```text
3fea20be042106fb5fe452ebe241dc5c3316eed6a893b00bf7ca2bcc0bef1b70
```

The continued two-move search is:

```bash
.venv/bin/python scripts/search_quartic_lr_two_move_circuits.py \
  --source-degree 2 \
  --source-parameter-bound 4 \
  --jet-max-exponent 12 \
  --workers 4 \
  --output artifacts/generated-results/quartic_lr_two_move_circuits.json
```

Its generated record has SHA-256

```text
28c9e1c2ed9c765fef7c51d7e8ace3262c6fac337c2b11d723f3dccdb3781826
```

The theorem-level affine and intrinsic checks are reproduced by

```bash
make verify-linear-torus-free
make verify-algebraic-torus-free
```
