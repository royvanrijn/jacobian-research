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

This statement holds in every degree `n>=3`. The proof is internal to the
polynomial-graph family and does not infer a uniform result from the finite
computational audit.

## Uniform genericity proof

Work first with degree-`n` polynomials modulo affine-linear summands. Adding
`aW+b` merely translates the dual coordinates by `(a,-b)`, so this coefficient
space has dimension `n-1` and retains all dual singularity data.

Every bad affine singularity forces one of three common-tangent contact
factorizations, with the displayed contact points distinct:

| Bad event | Factorization of `H(W)-ell(W)` |
|---|---|
| higher cusp at `r` | `(W-r)^4 Q(W)` |
| cusp at `r` plus another branch `u` | `(W-r)^3(W-u)^2 Q(W)` |
| three normalization points over one image | `(W-r)^2(W-u)^2(W-v)^2 Q(W)` |

If the total contact order is `m` at `k` marked points, then `Q` contributes
`n-m+1` coefficients and the points contribute `k` parameters. Each row
therefore has dimension

\[
k+(n-m+1)=n-2
\]

when it is nonempty; in smaller degrees it is empty. Compactifying the marked
points shows that the closure of each incidence image still has dimension at
most `n-2`. Their finite union is consequently a proper closed subset of the
`n-1` dimensional coefficient space. This simultaneously excludes multiple
roots of `H''`, cusp-branch collisions, tritangents, and coincident node
images. Two distinct regular branches can never be tangent to each other,
because their tangent directions are `(1,r)` and `(1,u)`, with determinant
`u-r`.

It remains to show that this open set meets the admissible weighted-seed
slice. Let `G` be a generic good polynomial, choose a generic tangent point
`alpha`, and let `beta!=alpha` be a simple residual intersection of that
tangent with the graph. Put `d=beta-alpha` and

\[
H(W)=G(\alpha+dW)-G(\alpha)-G'(\alpha)dW.
\]

Then

\[
H(0)=H'(0)=H(1)=0,
\qquad H'(1)\ne0.
\]

Moreover, writing `(s_G,t_G)=(G',WG'-G)` at `alpha+dW`,

\[
(s_H,t_H)=
\bigl(d(s_G-G'(\alpha)),\ t_G-\alpha s_G+G(\alpha)\bigr).
\]

Thus tangent-chord normalization is an affine source reparameterization
followed by an invertible affine target change and preserves every dual
singularity type. The tangent-chord incidence is irreducible and dominates
the unrestricted coefficient space. The remaining weighted condition
`H''(1)/(-H'(1))!=-2` removes only a proper closed subset: the admissible audit
family simplifies to

\[
H_n(W)=W^2+W^3+\cdots+W^{n-1}-(n-2)W^n
\]

and has `H''(1)/(-H'(1))=-4n/3`. Hence the good open set meets the admissible
slice for every `n>=3`, proving the theorem uniformly.

## Degree and birational normalization

The map `nu_H` is birational onto its image in characteristic zero. Indeed,
where `H''(r)!=0`,

\[
{dt\over ds}={rH''(r)\over H''(r)}=r,
\]

so the normalization coordinate is recovered rationally as the tangent slope
of the discriminant curve. After homogenization the three parameterizing
forms have degree `n` and no common factor. A generic line therefore pulls
back to a degree-`n` divisor, proving that the image curve has degree `n` and
normalization `P^1`.

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

nodes. As a regression certificate, the exact audit also checks the global
Tjurina number

\[
2(n-2)+{(n-2)(n-3)\over2},
\]

which is the minimal value contributed by precisely these ordinary cusps and
nodes and independently rules out higher collisions in the deterministic
examples through degree ten.

## Generic surjectivity

At a smooth discriminant point, the inverse polynomial has partition
`2+1+...`. At an ordinary cusp it has `3+1+...`, and at an ordinary node it
has `2+2+1+...`. Therefore:

- `n=3`: the cusp is a triple-root omitted value;
- `n=4`: the node is a double-double omitted value;
- `n>=5`: every discriminant value retains at least one simple root.

Hence a generic admissible weighted seed of every inverse degree `n>=5` is
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

The verifier checks the uniform incidence-dimension formulas and tangent-chord
normalization identities symbolically. It also retains deterministic rational
regression seeds in every inverse degree from three through ten, including
squarefree cusp polynomials, bitangent counts, the smooth infinity jet, global
singular-scheme length, and the multiplicity dictionary.

## Reference

For the surrounding projective-duality machinery, see Aliaksandr Yuran,
*Plucker Formulas for Plane Algebraic Curves with a Given Newton Polygon*,
arXiv:2204.04626. It treats inflections and bitangents as cusps and nodes of
the dual curve and proves generic nodal-cuspidal behavior under suitable
Newton-polygon hypotheses. The theorem above instead uses the elementary
contact-incidence dimension argument specialized to polynomial graphs.
