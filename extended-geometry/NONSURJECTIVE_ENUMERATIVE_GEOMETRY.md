# Enumerative skeleton of the nonsurjective locus

Fix an inverse degree `N>=3`.  The irreducible components of the
nonsurjective seed locus are

\[
 \mathcal C_{a,b}=\overline{\mathcal E_{2^a3^b}},
 \qquad 2a+3b=N.
\]

This note extracts the enumerative consequences that follow formally from
the contact-strata theorem and the explicit component normalizations.  It
also separates them from projective seed-closure degrees, for which a
boundary base-scheme calculation is still required.

## 1. Number and dimensions of the components

Let

\[
 c_N=\#\{(a,b)\in\mathbb Z_{\ge0}^2:2a+3b=N\}.
\]

Writing `N=6q+r`, with `0<=r<6`, gives the quasipolynomial

\[
 \boxed{
 c_N=
 \begin{cases}
 q,&r=1,\\
 q+1,&r=0,2,3,4,5.
 \end{cases}}
 \tag{1}
\]

Equivalently,

\[
 \boxed{\sum_{N\ge0}c_Nx^N={1\over(1-x^2)(1-x^3)}.}             \tag{2}
\]

Put `epsilon_N=0` for even `N` and `epsilon_N=1` for odd `N`.  The
components can be indexed consecutively by

\[
 b_j=\epsilon_N+2j,
 \qquad
 a_j={N-3b_j\over2},
 \qquad 0\le j<c_N.                                             \tag{3}
\]

Their dimensions and codimensions form consecutive staircases:

\[
 \boxed{
 \dim\mathcal C_j=\left\lfloor{N\over2}\right\rfloor-1-j,
 \qquad
 \operatorname{codim}_{\mathcal A_N}\mathcal C_j
 =\left\lfloor{N-3\over2}\right\rfloor+j.}                    \tag{4}
\]

Thus there is exactly one component in every dimension occurring in (4).
The simultaneous dimension--codimension enumerator is

\[
 \boxed{
 \sum_{N\ge3}\sum_{2a+3b=N}
 x^N u^{a+b-1}v^{a+2b-2}
 ={1\over uv^2}
 \left(
 {1\over(1-uvx^2)(1-uv^2x^3)}-1-uvx^2
 \right).}                                                     \tag{5}
\]

Specializing `u=v=1` recovers (2), after deleting degrees zero and two.

## 2. Degrees of the normalization maps

For every component, the quotient-coordinate normalization is

\[
 \widetilde{\mathcal C}_{a,b}
 =\{\Phi(Q^2R^3)=0\}\cap D(DA)
 \longrightarrow\mathcal C_{a,b}.
\]

It is finite and birational, so

\[
 \boxed{\deg(\widetilde{\mathcal C}_{a,b}/\mathcal C_{a,b})=1.} \tag{6}
\]

If instead all double-root and triple-root support points are labelled, the
ordered-root hypersurface maps generically to the quotient normalization with
degree

\[
 \boxed{a!b!.}                                                   \tag{7}
\]

These are intrinsic function-field degrees and do not depend on a
compactification.

## 3. Pairwise and multiple component intersections

Take two components whose indices differ by `k>0`.  In `(a,b)` notation they
are

\[
 (a,b)\quad\hbox{and}\quad(a-3k,b+2k).
\]

The common coarsening having the largest possible support is

\[
 \boxed{\nu_{a,b;k}=2^{a-3k}3^b6^k.}                            \tag{8}
\]

Indeed, the shared double and triple atoms remain separate, while each of
the `k` exchanges

\[
 2+2+2=3+3=6
\]

creates one sixfold block.  Any other common coarsening has no more support
points.  Consequently the unique top contact type in the set-theoretic
intersection has

\[
 \boxed{
 \dim\mathcal E_{\nu_{a,b;k}}=a+b-2k-1.}                        \tag{9}
\]

The intersection has codimension `2k` in the lower-`b` component and
codimension `k` in the higher-`b` component.  Its complete support is still
the common-coarsening union

\[
 \mathcal C_{a,b}\cap\mathcal C_{a-3k,b+2k}
 =\bigcup_{\substack{2^a3^b\preceq\nu\\
                     2^{a-3k}3^{b+2k}\preceq\nu}}
   \mathcal E_\nu.                                               \tag{10}
\]

More generally, the intersection of any collection of components depends
set-theoretically only on its two extreme indices.  To see this, decompose
each block of a common coarsening of the extremes.  The two atom allocations
in that block differ by an integral number of moves `3*2 <-> 2*3`.  The
allowable numbers of moves form an integer interval.  The Minkowski sum of
these intervals contains every integer between zero and the total `k`, so
every intermediate atomic type also refines the block decomposition.

The number of distinct index intervals of gap `k` in degree `N` has the
generating function

\[
 \boxed{
 \sum_{N\ge0}\#\{(i,i+k):0\le i<i+k<c_N\}x^N
 ={x^{6k}\over(1-x^2)(1-x^3)}.}                                \tag{11}
\]

Keeping the dimension of the top contact type gives the two-variable form

\[
 \boxed{
 \sum_{N\ge0}\sum_{k\ge1}\sum_i
 x^N t^{\dim\mathcal E_{\nu_i}+1}y^k
 ={1\over(1-tx^2)(1-tx^3)}
   {ytx^6\over1-ytx^6}.}                                       \tag{12}
\]

Equations (10)--(12) count supports and dimensions.  They do not assign
scheme-theoretic intersection multiplicities.  The separate
[weighted partition-lattice theorem](PARTITION_LATTICE_GEOMETRY.md) computes
the generic minimal-intersection transverse algebra as a tensor product of
`k[epsilon]/(epsilon^2)` blocks, and hence assigns transverse length `2^k`.
Consequently the scheme-weighted refinement of (12) is

\[
 \boxed{
 \sum_{N,k,i}2^k x^N t^{\dim\mathcal E_{\nu_i}+1}y^k
 ={1\over(1-tx^2)(1-tx^3)}
   {2ytx^6\over1-2ytx^6}.}                                     \tag{12a}
\]

More finely, replacing `2^k` by the tangent-cone Hilbert polynomial
`(1+z)^k` replaces the last factor by

\[
 {y(1+z)tx^6\over1-y(1+z)tx^6}.                                \tag{12b}
\]

## 4. Normalization branches along the top intersections

Consider the interval of components from index `i` to index `i+k`, and its
top contact type (8).  For the component at offset `s`, with `0<=s<=k`, each
of the `k` distinct sixfold roots must be allocated either as three roots of
`Q` or as two roots of `R`.  Exactly `s` of them receive the latter
allocation.  Therefore

\[
 \boxed{
 r_{k,s}=\binom{k}{s},
 \qquad
 \sum_{s=0}^k r_{k,s}z^s=(1+z)^k.}                              \tag{13}
\]

This is the degree of the normalization fiber over the generic point of the
top contact stratum, not the generic degree (6) of the normalization map.
The endpoint components have one branch, while every interior component has
more than one branch when `k>=2`.  Hence

\[
 \boxed{
 \mathcal E_{\nu_{a,b;k}}
 \subseteq\operatorname{Sing}(\mathcal C_{a-3s,b+2s})
 \quad(0<s<k).}                                                  \tag{14}
\]

The first instance is `k=2`, `N=12`: the middle `(a,b)=(3,2)` component has
two branches over the `(6,6)` stratum.  Formula (13) recovers the branch
combinatorics behind the degree-twelve local calculation and predicts the
binomial branch profiles in higher multiple intersections.

For an arbitrary collision partition `nu=(m_1,...,m_e)`, the complete branch
enumerator remains

\[
 \boxed{
 r_{a,b}(\nu)=[U^aV^b]
 \prod_{\rho=1}^e
 \left(\sum_{2i+3j=m_\rho}U^iV^j\right).}                       \tag{15}
\]

Thus the collision partition stratifies every normalization fiber.  The
exact maximal stratum is smooth, so

\[
 \boxed{
 \bigcup_{\substack{2^a3^b\preceq\nu\\r_{a,b}(\nu)>1}}
 \mathcal E_\nu
 \subseteq\operatorname{Sing}(\mathcal C_{a,b})
 \subseteq\mathcal C_{a,b}\setminus\mathcal E_{2^a3^b}.}      \tag{16}
\]

The left inclusion is generally strict: one-branch ramification and
non-normal cusp behavior require differential or conductor calculations.

## 5. The projective coincident-root Chow package

There is one compact enumerative package available before projection to
normalized seed coefficients.  Let

\[
 \mu_{a,b}:\mathbb P^a\times\mathbb P^b\longrightarrow\mathbb P^N,
 \qquad(Q,R)\longmapsto Q^2R^3,
\]

and let `h` be the hyperplane class on `P^N`, while `u,v` are the two source
hyperplane classes.  Unique factorization makes `mu_(a,b)` generically
degree one, and

\[
 \mu_{a,b}^*h=2u+3v.                                            \tag{17}
\]

For the projective coincident-root locus `Delta_(a,b)` this gives

\[
 \boxed{
 \deg\Delta_{a,b}
 =\binom{a+b}{a}2^a3^b,
 \qquad
 [\Delta_{a,b}]
 =\binom{a+b}{a}2^a3^b\,h^{N-a-b}.}                             \tag{18}
\]

The degree generating function is especially simple:

\[
 \boxed{
 \sum_{a,b\ge0}\deg(\Delta_{a,b})
 x^{2a+3b}t^{a+b}
 ={1\over1-2tx^2-3tx^3}.}                                      \tag{19}
\]

Pulling back one projective linear incidence condition gives the total
divisor class `2u+3v`.  Its source multidegrees are

\[
 \boxed{
 \int_{\mathbb P^a\times\mathbb P^b}
 u^pv^q(2u+3v)^{r+1}
 =\binom{r+1}{a-p}2^{a-p}3^{b-q},}                              \tag{20}
\]

whenever `p+q+r=a+b-1`, and are zero outside the natural ranges.  The total
hyperplane-section Chow classes have the formal series

\[
 \sum_{a,b\ge0}
 \binom{a+b}{a}2^a3^b
 x^{2a+3b}z^{a+2b+1}
 ={z\over1-2zx^2-3z^2x^3}.                                     \tag{21}
\]

These formulas are useful inputs, but (18)--(21) are **not** yet degrees or
Chow classes of the closures inside the normalized seed space.

## 6. Why the seed-closure degree is a separate problem

The normalized seed remembers the coefficients of `M=Q^2R^3` only modulo
the affine span of `1,W`, followed by division by

\[
 D=M'(1)-M'(0).
\]

On the monic affine chart this is finite and gives the normalization theorem.
On the binary-form compactification it is a rational linear projection with
a boundary base scheme.  Therefore its image degree is the moving
intersection number after subtracting the Segre contribution of that base
scheme, not the raw number (18).

Degree four already forces this distinction.  The binary square locus has
degree `2^2=4`, whereas the nonsurjective normalized quartic locus is dense
in the one-dimensional seed space, so its projective seed closure has degree
one.

For fixed `(a,b)`, the exact next calculation is therefore:

1. compactify `V(Phi)` in `P^a x P^b` and saturate away boundary components;
2. form the Rees algebra of the retained-coefficient base ideal;
3. compute its Segre class, or equivalently the multidegree of the graph;
4. push the graph class to the projective normalized seed hyperplane.

This produces the actual degrees, multidegrees, and Chow classes requested
for component closures.  The raw rational series (19)--(21) is the
uncorrected term; the missing series is precisely the boundary Segre-class
correction.  Minimal pairwise intersection weights are computed in the
partition-lattice theorem; nonminimal collision blocks and complete
multi-sheet conductor equalizers still require their normal cones and
conductors.

## Executable certificate

Run

```bash
.venv/bin/python scripts/verify_nonsurjective_enumerative_geometry.py
```

The script checks the quasipolynomial and all rational generating functions,
the component staircases, the interval-intersection theorem through degree
twenty-two, the binomial normalization fibers, and the coincident-root Chow and
multidegree formulas.
