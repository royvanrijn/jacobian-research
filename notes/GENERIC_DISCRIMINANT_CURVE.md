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

Fix an algebraically closed field `k` of characteristic zero.  Work in the
affine coefficient space

\[
\mathcal B_n=\operatorname{Spec}k[h_2,\ldots,h_n]\simeq\mathbb A^{n-1},
\qquad H(W)=\sum_{j=2}^n h_jW^j.
\]

This is the space of degree-at-most-`n` polynomials modulo affine-linear
summands. Adding `aW+b` translates the dual coordinates by `(a,-b)` and does
not change any singularity type. The exact-degree locus is the open set
`B_n^o=D(h_n)`.

### Projective normalization and its point at infinity

For `[R:U] in P^1`, the projective extension of the tangent-line map is

\[
\overline\nu_H[R:U]=[S_H(R,U):T_H(R,U):Z(R,U)],
\]

where the three degree-`n` forms are

\[
\begin{aligned}
S_H&=\sum_{j=2}^n jh_jR^{j-1}U^{n-j+1},\\
T_H&=\sum_{j=2}^n (j-1)h_jR^jU^{n-j},\\
Z&=U^n.
\end{aligned}
\]

On `U=1` this is `[H'(R):RH'(R)-H(R):1]`.  On `B_n^o` the forms have no
common zero.  The only parameter point over the target line at infinity is
`[R:U]=[1:0]`, and its image is `[0:1:0]`.  In the local coordinate `q=U/R`
and the target chart `T=1`,

\[
{S_H\over T_H}={n\over n-1}q+O(q^2),\qquad
{Z\over T_H}={1\over(n-1)h_n}q^n+O(q^{n+1}).
\]

Thus `q` is recovered from the first target coordinate.  Uniformly on
`D(h_n)`, the map is a closed immersion in a formal neighborhood of infinity.
In particular the point at infinity is smooth, cannot meet a finite branch,
and cannot be the limit of two distinct normalization points with a common
image.

### The three compactified bad-incidence varieties

For a contact pattern

\[
\mu=(\mu_1,\ldots,\mu_k),\qquad
m=\sum_i\mu_i,\qquad d=n-m,
\]

put

\[
M_\mu(W)=\prod_{i=1}^k(W-r_i)^{\mu_i}.
\]

When `d>=0`, an affine common-tangent incidence has the coefficient identity

\[
H(W)-\ell(W)=M_\mu(W)Q(W),\qquad \deg Q\le d. \tag{1}
\]

Since `ell` is affine-linear, equation (1) says exactly that the coefficients
of `W^2,...,W^n` in `H` and `M_mu Q` agree.  Consequently it defines a
polynomial map

\[
\theta_\mu:\mathbb A^k\times\mathbb A^{d+1}\longrightarrow\mathcal B_n.
\tag{2}
\]

The exact-contact incidence is the open subset where the marked roots are
distinct, `Q(r_i)!=0`, and the product has degree `n`.

Here is an explicit proper compactification of (2).  Write a marked root as
`[R_i:U_i]`, put

\[
L_i(W,V)=U_iW-R_iV,\quad
D=\prod_iU_i^{\mu_i},\quad
M_\mu^h=\prod_iL_i^{\mu_i},
\]

and write

\[
Q^h(W,V)=\sum_{a=0}^d q_aW^aV^{d-a}.
\]

Compactify the residual coefficients as
`[q_infty:q_0:...:q_d] in P^(d+1)` and the target as
`[z:h_2:...:h_n] in P^(n-1)`.  If `P_j` denotes the coefficient of
`W^jV^(n-j)` in `M_mu^h Q^h`, define

\[
\overline{\mathcal I}_\mu\subset
(\mathbb P^1)^k\times\mathbb P^{d+1}\times\mathbb P^{n-1} \tag{3}
\]

to be the graph closure, equivalently the subscheme obtained from

\[
Dq_\infty h_j-zP_j=0\qquad(2\le j\le n) \tag{4}
\]

by saturation with respect to `z D q_infty`.  This saturation is part of the
definition: it removes components created solely by homogenizing away from
the affine graph chart.  If `d<0`, the incidence is empty.

The three bad patterns and their compactifications are

| Pattern `mu` | Event | `m` | `k` | `d+1` residual coefficients |
|---|---|---:|---:|---:|
| `(4)` | nonordinary ramified branch | 4 | 1 | `n-3` |
| `(3,2)` | ramified branch sharing its image | 5 | 2 | `n-4` |
| `(2,2,2)` | at least three normalization points over one image | 6 | 3 | `n-5` |

The source of every nonempty affine graph in (2) is irreducible of dimension

\[
k+(d+1)=k+n-m+1=n-2. \tag{5}
\]

Its graph closure (3) is therefore irreducible of the same dimension; taking
a closure introduces no new irreducible component.  Projection from (3) to
`P^(n-1)` is proper, so its image is closed and is exactly the Zariski closure
of the corresponding affine incidence image.  Hence

\[
\dim\overline{\theta_\mu(\mathcal I_\mu^o)}\le n-2. \tag{6}
\]

This proves the dimension bound for the closure, not merely for the
distinct-root chart.  Root collisions, residual factors meeting marked
roots, and every projective boundary stratum remain inside the same
irreducible graph closure and cannot raise its dimension.

### Exhaustion of every bad affine singularity

For a finite parameter `r`, let `ell_r` be the tangent line to the graph of
`H` at `r`.  Its contact order is

\[
c(r)=\operatorname{ord}_{W=r}(H(W)-\ell_r(W))\ge2.
\]

The derivative formula

\[
\nu_H'(r)=H''(r)(1,r) \tag{7}
\]

gives the complete branch classification:

1. `c(r)=2` exactly when the normalization is immersive at `r`.
2. `c(r)=3` exactly when `H''` has a simple zero at `r`.  Then
   `det(nu_H''(r),nu_H'''(r))=2H'''(r)^2!=0`, so the image branch is an
   ordinary cusp.
3. `c(r)>=4` gives pattern `(4)` after absorbing all excess contact into `Q`.

Now suppose a singular image point has more than one normalization preimage.
Equality of two images says precisely that one affine line is tangent at both
graph points, so every preimage contributes contact at least two to (1).
There are only three possibilities not already ordinary:

- if one of exactly two branches is ramified, their contact orders dominate
  `(3,2)`; this also includes cusp-cusp collisions and higher contact after
  moving the excess multiplicity into `Q`;
- if there are at least three preimages, any three give `(2,2,2)`; this also
  includes coincident node images and fibers with four or more preimages;
- if there are exactly two unramified preimages `r!=u`, the two branch tangent
  directions are `(1,r)` and `(1,u)`, whose determinant is `u-r!=0`.  The
  singularity is therefore an ordinary transverse node, never a tacnode.

Together with the one-branch classification, these alternatives prove that
`(4)`, `(3,2)`, and `(2,2,2)` exhaust every affine singularity worse than an
ordinary node or ordinary cusp.  There is no additional nonbirational case:
on `H''!=0`, the function `dt/ds=r` recovers the normalization parameter from
the image function field.

### Complete boundary analysis

For clarity, every degeneration of the open incidences is listed here.

1. **Finite marked-root collisions.**  Merging marked roots adds their contact
   orders: `(3,2)` can specialize to contact at least five, while `(2,2,2)`
   can specialize to `(4,2)` or contact at least six.  These are points of the
   same closures and are already dominated by one of the three patterns.
2. **Residual collision.**  If `Q(r_i)=0`, the corresponding contact order
   increases.  This is again excess multiplicity inside the same closure; it
   includes patterns such as `(3,3)`, `(4,2)`, and higher cusps.
3. **A marked point tends to parameter infinity.**  A ramified marked point
   would make the homogenization of `H''` vanish at `[1:0]`, but its value
   there is `n(n-1)h_n`; hence this forces degree drop `h_n=0`.  An unramified
   point at infinity cannot share an image with a finite point because their
   target `Z`-coordinates differ.  Two or more marked points cannot collide
   there off the diagonal because the local expansion above makes
   `overline(nu)_H` a closed immersion at infinity.
4. **Residual coefficients tend to infinity.**  Once all marked roots remain
   finite, a nonzero limiting direction `Q_infty` would have to make every
   coefficient of degree at least two in `M_mu Q_infty` vanish in order for
   the target coefficients to stay affine and bounded.  Thus
   `M_mu Q_infty` would be affine-linear, contradicting
   `deg M_mu=m>=4`.
5. **The coefficient target tends to infinity.**  This is exactly the boundary
   `z=0` of `P^(n-1)` and does not lie over `B_n`.
6. **Degree drop.**  The remaining boundary is contained in `h_n=0`, outside
   `B_n^o`.  After cancelling the common power of `U`, it is the analogous
   lower-degree problem and supplies no singularity of an exact degree-`n`
   discriminant.

These cases exhaust the boundary of (3): it is the union of `z=0`,
`q_infty=0`, some `U_i=0`, the finite diagonals, the residual-resultant
divisors, and `h_n=0`.  Thus no unexamined infinity component projects into
the exact-degree affine coefficient space.

### The good open set

Let `Z_4`, `Z_32`, and `Z_222` be the three closed projective images from
(3).  By (6), each has dimension at most `n-2`, whereas `B_n` has dimension
`n-1`.  Their union with the degree-drop divisor is a proper closed subset.
Its complement `U_n` has squarefree `H''`, no cusp-branch collision, and no
tritangent or coincident node images.  The exhaustion proof and the infinity
analysis show that every discriminant belonging to `U_n` has only ordinary
cusps and ordinary nodes and is smooth at infinity.

### Why the admissible weighted slice meets the good open

It remains to prove that `U_n` is not lost after imposing the weighted endpoint
conditions.  Define the tangent-chord incidence

\[
\mathcal T_n=\left\{(G,\alpha,\beta):\alpha\ne\beta,\quad
{G(\beta)-G(\alpha)-G'(\alpha)(\beta-\alpha)
 \over(\beta-\alpha)^2}=0\right\}. \tag{8}
\]

The divided expression is polynomial.  Its coefficient of `g_2` is exactly
one, so (8) solves uniquely for `g_2`.  Consequently

\[
\mathcal T_n\simeq
(\mathbb A^2\setminus\Delta)\times\mathbb A^{n-2} \tag{9}
\]

and is irreducible.  The projection `T_n -> B_n` is dominant: for an
exact-degree `G` and a generic `alpha`, the divided tangent-intersection
polynomial has degree `n-2` in `beta`, and it does not vanish at `beta=alpha`
because that value is `G''(alpha)/2`.  Over the algebraic closure it therefore
has a distinct residual root.

Put `d=beta-alpha` and

\[
H(W)=G(\alpha+dW)-G(\alpha)-G'(\alpha)dW. \tag{10}
\]

Then `H(0)=H'(0)=H(1)=0`.  Writing `(s_G,t_G)=(G',WG'-G)` at
`alpha+dW`, one obtains

\[
(s_H,t_H)=
\bigl(d(s_G-G'(\alpha)),\ t_G-\alpha s_G+G(\alpha)\bigr). \tag{11}
\]

Equation (11) is an affine reparameterization of `P^1` followed by an
invertible affine target transformation of determinant `d`.  It preserves
the projective discriminant, including the local type at infinity.  Hence the
inverse image of `U_n` in `T_n` is a nonempty open subset.

The remaining weighted conditions are also a nonempty open subset of `T_n`:

\[
\deg H=n,\qquad H'(1)\ne0,\qquad H''(1)-2H'(1)\ne0. \tag{12}
\]

Nonemptiness is witnessed on (8) by `G=H_n`, `alpha=0`, `beta=1`, where

\[
H_n(W)=W^2+W^3+\cdots+W^{n-1}-(n-2)W^n
\]

and

\[
{H_n''(1)\over-H_n'(1)}=-{4n\over3}\ne-2.
\]

Because `T_n` is irreducible, its nonempty good open and its nonempty
weighted-admissible open intersect.  Equation (10) at a point of this
intersection is a good admissible seed.  Finally, scaling by the nonzero
constant `-1/H'(1)` reaches the normalized slice `H'(1)=-1`; the ratio in
(12) and every dual singularity type are unchanged.  Therefore the normalized
admissible weighted slice meets `U_n` for every `n>=3`.

## Universal saturated incidence ideals

The module `jcsearch.discriminant_geometry` now represents the theorem's bad
loci directly over the universal admissible chart

\[
H(W)=\sum_{k=3}^n h_k(W^k-W^2).
\]

This imposes `H(0)=H'(0)=H(1)=0` identically. The coefficient open set is
represented by the single factor

\[
\Omega_n=h_nH'(1)\bigl(H''(1)-2H'(1)\bigr),
\]

excluding degree drop, `c=-H'(1)=0`, and the forbidden weighted value
`H''(1)/c=-2`.

Write `D_s(r,u),D_t(r,u)` for the two diagonal-divided bitangent equations.
The reusable incidence ideals and saturation factors are:

| Incidence | Equations | Saturation factor besides `Omega_n` |
|---|---|---|
| ordinary ordered bitangent | `D_s(r,u),D_t(r,u)` | `(r-u)H''(r)H''(u)` |
| higher cusp | `H''(r),H'''(r)` | `1` |
| ordinary cusp plus branch | `H''(r),D_s(r,u),D_t(r,u)` | `(r-u)H'''(r)H''(u)` |
| tritangent | `D(r,u),D(r,v)` | `(r-u)(r-v)(u-v)H''(r)H''(u)H''(v)` |

If `S` is the full displayed saturation factor, the implementation adjoins a
Rabinowitsch variable `z` and the equation `1-zS=0`. Eliminating the marked
points and `z` therefore produces the closed coefficient-space stratum without
retaining diagonals, nonordinary cusp endpoints, or inadmissible parameters.
`incidence_elimination_generators` performs this elimination when a concrete
degree is small enough for a direct Groebner calculation.

The executable audit verifies every contact factorization and saturation gate
symbolically. As a genuine elimination regression, the degree-four
higher-cusp ideal recovers exactly the discriminant of `H''` in the admissible
coefficient space.

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
.venv/bin/python scripts/verify_universal_discriminant_incidences.py
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
