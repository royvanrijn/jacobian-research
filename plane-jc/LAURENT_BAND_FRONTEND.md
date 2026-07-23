# Laurent-polygon to de Rham front end

## Status

[`cas/laurent_band_frontend.py`](cas/laurent_band_frontend.py) now implements
the exact middle front end

```text
published Laurent polygons
    -> lattice points
    -> unimodular monomial chart
    -> z-band supports
    -> bracket layers
    -> WeightedWronskianIR
    -> de Rham compiler
```

It deliberately does not implement

```text
admissible corner chain -> published Laurent polygons.
```

That arrow uses the normal-form, predecessor/successor, approximate-root, and
Laurent-translation theorems of the plane-JC reduction.  It is not determined
by the corner list alone.

The code represents this boundary explicitly as a
`NewtonNormalFormCertificate`.  A certificate maps one chain to zero or more
Laurent cases, records the transformations and theorem source, and separately
records whether the branch list is exhaustive.  Compilation is permitted only
when the list is nonempty, exhaustive, and has no missing data.

## 1. Exact support compilation

For the audited chart

\[
 t=xy^2,\qquad z=y^{-1},
\]

the inverse monomial data are

\[
 x=tz^2,\qquad y=z^{-1},
\]

so the exponent map is

\[
 (i,j)\longmapsto(i,2i-j).
\]

The exponent matrix is unimodular and `[t,z]_(x,y)=-1`.  The front end
enumerates every lattice point in each convex polygon and groups the images
by their `z` exponent.  It reproduces:

| Case | `P` lattice points | `P` bands | `Q` lattice points | `Q` bands |
| --- | ---: | --- | ---: | --- |
| 1 | 61 | `2,1,0,...,-8` | 125 | `3,2,1,...,-12` |
| 2 | 25 | `2,1,0` | 47 | `3,2,1,0` |

Both top supports are

\[
 A:\ t^1,\ldots,t^8,
 \qquad
 D:\ t^2,\ldots,t^{12}.
\]

The normalization fixes the coefficients of `t` and `t^8` in `A` to one.
Since `x^2=t^2z^4`, the top bracket layer compiles automatically to

\[
 2AD'-3A'D=t^2.
\]

Feeding that result to the de Rham compiler gives the certified genus-three,
six-coordinate obstruction for either polygon.

## 2. Bracket-layer compiler

For constant `J=[t,z]_(x,y)`, the implementation uses

\[
 [P_i(t)z^i,Q_j(t)z^j]_{x,y}
 =Jz^{i+j-1}\bigl(jP_i'Q_j-iP_iQ_j'\bigr).
\]

On the audited upper bands it reproduces all five equations `J4,...,J0`,
including their signs.  Zero formal layers such as the bracket of the two
`z^0` bands are removed.  The current implementation intentionally requires
a unimodular Laurent chart and constant `J`; a nonconstant Jacobian monomial
would also shift band degrees and needs a larger IR.

## 3. Why a corner chain is insufficient

The audited example itself proves that the map from chains to polygon pairs
is not a function of the printed corners alone.  The single chain

\[
 (8,28)\longrightarrow(11/4,7)
\]

produces two distinct Proposition 4.3 polygon pairs.  Their common upper
bands happen to compile to the same first block, while their lower bands lead
to different later systems.

Therefore a compiler may not infer a unique Laurent polygon from a row of the
2017 complete-chain table.  It needs a normal-form certificate containing at
least:

- the sequence of flips, Laurent translations, and monomial maps;
- every polygon branch and the theorem making the branch list exhaustive;
- exact vertex nonvanishing and coefficient normalizations;
- the full lattice support or inequalities defining it;
- the transformed bracket monomial.

Once those items are supplied, the new front end performs the remaining
support and bracket transcription mechanically.

There is now a complementary
[log-boundary compiler](LOG_BOUNDARY_COMPILER.md).  It consumes certified
positive local branch scales and builds the toroidal proximity/boundary
package.  It does not weaken the refusal above: neither polygons nor positive
local scales are functions of the printed corners alone.  A normal-form
theorem must supply both exhaustive branch lists before the Laurent-band and
log-boundary compilers can be used together.

When a compiled Poisson-square truncation lands in the three-layer degree
box `(3,4;2,3)`, the
[Poisson-square rigidity theorem](POISSON_SQUARE_RIGIDITY.md) replaces its
geometric reduced coefficient locus by three classified components: the
four-parameter tangent closure, `C=0`, and `A=0`.  Only the nonreduced
intersection structure remains a separate input; the three generic
multiplicities are already `2,3,1`.

This is also an executable distinction.  The audited `(72,108)` certificate
contains exactly two distinct cases and compiles both.  The `(75,125)` F2
record contains no cases, is marked non-exhaustive, and raises an error if a
caller attempts to compile it.

## 4. Consequence for `(75,125)`

Family `F2` at `j=1` supplies only

\[
 (5,20)\longrightarrow(7/5,2),\qquad(m,n)=(3,5).
\]

The chain theorem nevertheless forces a useful partial certificate.  With
`X=x^(1/5)`, the required Puiseux translation is
`y -> y+lambda/X`, the bracket becomes `X^4`, and two consecutive nonzero
edges are fixed.  The terminal type-I edge normalizes uniquely and has no
residual de Rham obstruction.  These exact facts and the remaining proof
obligations are derived in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md) and encoded by
[`cas/f2_75_125_frontend.py`](cas/f2_75_125_frontend.py).

The older `(50,75)` calculation concerns the `j=0`, `(2,3)` member and begins
with two additional cases `gamma=2,3`; it explicitly says that the first
normal-form part is not proved there because it is used only to revisit a
known case.  Consequently it is neither a proof nor a template for the
`(3,5)` member.

The next genuine frontier task is therefore the lower-boundary part of a new
F2 `j=1` normal-form proposition: prove the `gamma` branches exhaustive and
produce the complete Laurent polygon pairs.  Until that theorem is derived,
`(75,125)` has no honest complete band input for the de Rham engine.

## 5. Reproduction

Run:

```bash
python3 plane-jc/cas/test_laurent_band_frontend.py
```

Expected output:

```text
PASS: both audited polygons compile to their exact Laurent band supports
PASS: one chain expands to two exhaustive, distinct Laurent cases
PASS: polygon top bands compile to the certified genus-three first block
PASS: the front end reproduces all five upper bracket layers J4,...,J0
PASS: an unsupported F2 frontier cannot compile as a normal form
```
