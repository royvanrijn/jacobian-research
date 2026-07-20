# Generic discriminant curve as the dual of a polynomial graph

For an inverse primitive `H` of degree `n`, the repeated-root normalization is

\[
\nu_H(r)=\bigl(H'(r),rH'(r)-H(r)\bigr).
\]

This is precisely the parameter space of tangent lines

\[
Y=H(r)+H'(r)(W-r)=sW-t
\]

to the graph `Y=H(W)`. Thus the discriminant is the projective dual curve of
the polynomial graph, and its singularities record exceptional tangent lines.

## Generic nodal-cuspidal theorem

On a Zariski-open subset of admissible degree-`n` primitives:

- the projective discriminant is a rational curve of degree `n`;
- its unique point at infinity is smooth;
- it has exactly `n-2` ordinary cusps;
- it has exactly
  \[
  { (n-2)(n-3)\over2}
  \]
  ordinary nodes;
- it has no other singularities.

The excluded parameter conditions are algebraic: a multiple root of `H''`, a
cusp sharing its tangent-line image with another normalization point, a
tritangent, coincident node images, nontransverse self-intersection, or a
degeneration of the infinity jet. They are obtained by eliminating the
normalization variables from the corresponding incidence equations after
saturating away diagonals and cusp factors. Thus each condition is Zariski
closed.

The classical generic-duality theorem says that the dual of a generic plane
curve has only ordinary nodes and cusps. Applying its incidence/transversality
argument to this polynomial-graph family proves that the complement of the bad
loci above is open. Nonemptiness is certified computationally here, over the
rationals, in each degree `3,...,10`. Extending this implementation-level
nonemptiness certificate uniformly to every degree still requires either a
uniform explicit seed family proof or a graph-family-specific transversality
proof; the verifier does not silently infer all degrees from ten examples.

## Cusps are flexes

Differentiating gives

\[
\nu_H'(r)=H''(r)(1,r).
\]

Hence ramification occurs exactly at the inflection points `H''(r)=0`. If the
root is simple, then at that point

\[
\nu_H''(r)=H'''(r)(1,r)
\]

and

\[
\det\bigl(\nu_H''(r),\nu_H'''(r)\bigr)=2H'''(r)^2\ne0.
\]

It is therefore an ordinary cusp. A generic `H''` is squarefree of degree
`n-2`, giving exactly `n-2` cusps.

## Nodes are bitangents

Two distinct normalization points `r,u` have the same image precisely when
the graph has one tangent line at both points. After removing the diagonal,
the exact bitangent equations are

\[
{H'(r)-H'(u)\over r-u}=0,
\]

\[
{rH'(r)-H(r)-uH'(u)+H(u)\over r-u}=0.
\]

They must be saturated by

\[
(r-u)H''(r)H''(u)
\]

to exclude the diagonal and cusp endpoints. At a genuine bitangent, the two
branch directions are `(1,r)` and `(1,u)`. Their determinant is `u-r`, so the
intersection is automatically transverse and produces an ordinary node.

The implementation passes to the symmetric coordinates `r+u,ru`, so each
unordered bitangent is counted once. Exact resultant calculations give

\[
{(n-2)(n-3)\over2}
\]

solutions for the deterministic audit family.

The resultant retains `n-2` diagonal flex solutions before removal. Its
remaining factor has degree `(n-2)(n-3)/2` and is squarefree and coprime to the
cusp-image factor. Restoring the two orderings gives saturated ordered-pair
quotient length `(n-2)(n-3)`. The global singular-scheme calculation below
then rules out hidden triple fibers and coincident node images in each audited
example.

## The point at infinity

Homogenize the parameterization as

\[
[H'(r):rH'(r)-H(r):1].
\]

There is one point at infinity, `[0:1:0]`. If `q=1/r` and `h_n` is the leading
coefficient of `H`, then in the chart where the middle coordinate is one,

\[
{H'(r)\over rH'(r)-H(r)}={n\over n-1}q+O(q^2),
\]

\[
{1\over rH'(r)-H(r)}={1\over(n-1)h_n}q^n+O(q^{n+1}).
\]

The first coordinate is a uniformizer, proving smoothness at infinity.

## Genus count

A plane curve of degree `n` has arithmetic genus

\[
p_a={(n-1)(n-2)\over2}.
\]

The normalization is rational, so its geometric genus is zero. Each ordinary
cusp and node has delta invariant one. Subtracting the `n-2` cusp
contributions leaves

\[
{(n-1)(n-2)\over2}-(n-2)
={(n-2)(n-3)\over2}
\]

nodes. The exact audit also checks the global Tjurina number

\[
2(n-2)+{(n-2)(n-3)\over2},
\]

which is the minimal value contributed by precisely these ordinary cusps and
nodes and rules out higher collisions in the deterministic examples.

## Generic surjectivity

At a smooth discriminant point, the inverse polynomial has partition
`2+1+...`. At an ordinary cusp it has `3+1+...`, and at an ordinary node it
has `2+2+1+...`. Therefore:

- `n=3`: the cusp is a triple-root omitted value;
- `n=4`: the node is a double-double omitted value;
- `n>=5`: every discriminant value retains at least one simple root.

Hence a generic admissible weighted seed of inverse degree `n>=5` is
surjective over the algebraic closure.

For inverse degree five, the exceptional equation `F(R,P)=0` from
`OMITTED_VALUE_CLASSIFICATION.md` now has a geometric meaning: it is the locus
where a cusp image acquires another normalization preimage. The inverse
partition becomes `3+2`, exhausting the degree and creating an omitted value.

## Multiplicity dictionary

| Root partition | Dual-curve geometry |
|---|---|
| `3,1,...` | ordinary cusp |
| `2,2,1,...` | ordinary node |
| `4,1,...` | higher cusp |
| `3,2,1,...` | cusp branch meeting another tangent branch |
| `2,2,2,1,...` | tritangent / triple normalization point |
| `n` | maximally ramified cusp |

The omitted-value classifier can therefore be read geometrically: enumerate
the discriminant singularities whose normalization-contact multiplicities
exhaust the full inverse degree.

Run:

```bash
.venv/bin/python scripts/verify_generic_discriminant_geometry.py
```

The verifier checks deterministic rational admissible seeds in every inverse
degree from three through ten, including squarefree cusp polynomials,
bitangent counts, the smooth infinity jet, global singular-scheme length, and
the multiplicity dictionary.

## Reference

For the surrounding projective-duality machinery, see Aliaksandr Yuran,
*Plucker Formulas for Plane Algebraic Curves with a Given Newton Polygon*,
arXiv:2204.04626. It treats inflections and bitangents as cusps and nodes of
the dual curve and proves generic nodal-cuspidal behavior under suitable
Newton-polygon hypotheses. The exact graph-family checks in this repository
are deliberately stated separately from those hypotheses.
