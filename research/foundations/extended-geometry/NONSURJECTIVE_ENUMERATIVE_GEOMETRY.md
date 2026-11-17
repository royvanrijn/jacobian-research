# Enumerative skeleton of the nonsurjective locus

Fix an inverse degree `N>=3`.  The irreducible components of the
nonsurjective seed locus are

\[
 \mathcal C_{a,b}=\overline{\mathcal E_{2^a3^b}},
 \qquad 2a+3b=N.
\]

This note extracts the enumerative consequences that follow from the
contact-strata theorem and the explicit component normalizations.  The
retained-coefficient projection has one boundary base point; a weighted
local-intersection calculation below gives its exact Segre correction and
hence the projective degrees and full graph Chow classes of all normalized
seed-component closures.

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
scheme-theoretic intersection multiplicities.  The scheme-weighted and
tangent-cone refinements of (12) lie outside the scope of this enumerative
skeleton; they are stated and proved as a corollary of the separate
[weighted partition-lattice theorem](PARTITION_LATTICE_GEOMETRY.md#weighted-interval-series).

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

These formulas are the uncorrected inputs for the retained-coefficient
projection.  They are not by themselves degrees or Chow classes of the
closures inside the normalized seed space.

## 6. Exact projective degrees after retained-coefficient projection

The normalized seed remembers the coefficients of `M=Q^2R^3` only modulo
the affine span of `1,W`, followed by division by

\[
 D=M'(1)-M'(0).
\]

On the monic affine chart this is finite and gives the normalization theorem.
On the binary-form compactification it is a rational linear projection with
one boundary base point.  The image degree is therefore the moving
intersection number after subtracting the local Segre contribution at that
point, not the raw number (18).

Degree four already forces this distinction.  The binary square locus has
degree `2^2=4`, whereas the nonsurjective normalized quartic locus is dense
in the one-dimensional seed space, so its projective seed closure has degree
one.

Put `ell=a+b`.  The exact answer is

\[
 \boxed{
 \deg\mathcal C_{a,b}
 =\binom{\ell}{a}\left(2^a3^b-(\ell+1)\right).}                \tag{22}
\]

Equivalently,

\[
 \boxed{
 \sum_{a,b\ge0}\deg\mathcal C_{a,b}\,
 t^{a+b}x^{2a+3b}
 ={1\over1-2tx^2-3tx^3}
  -{1\over(1-tx^2-tx^3)^2}.}                                  \tag{23}
\]

The terms `(a,b)=(0,0),(1,0)` on the right cancel to zero, so (23) may be
written over all nonnegative pairs even though the exceptional components
start in inverse degree three.

### Projective compactification and its unique base point

Write binary forms in variables `(W,V)` and set

\[
 S=\mathbb P^a\times\mathbb P^b,
 \qquad L=\mathcal O_S(2,3).
\]

The multiplication morphism

\[
 \mu:S\longrightarrow\mathbb P(k[W,V]_N),
 \qquad([Q],[R])\longmapsto[Q^2R^3]
\]

satisfies `mu^*O(1)=L`.  The linear functional

\[
 \Phi(M)=\sum_{j=2}^N m_j
\]

vanishes on the projection center

\[
 \mathbb P\langle V^N,WV^{N-1}\rangle.
\]

Let `X=mu^*V(Phi)`.  The retained coefficients `m_2,...,m_N` define on `X`
a rational map to the projective normalized seed space

\[
 X\dashrightarrow
 \mathbb P\!\left(\ker(\Phi)/
 \langle V^N,WV^{N-1}\rangle\right)
 \simeq\mathbb P^{N-3}.                                      \tag{24}
\]

Its set-theoretic base locus consists of the single point

\[
 p=([V^a],[V^b]).                                             \tag{25}
\]

Indeed, every nonconstant member of
`<V^N,WV^(N-1)>` has a simple root, whereas every root multiplicity of
`Q^2R^3` has the form `2i+3j` and is therefore either zero or at least two.
The remaining center member is `V^N`, whose factorization gives exactly
(25).

For `ell>=2`, `X` has no component in either leading-coefficient boundary
hyperplane.  For example, on a double-factor boundary one may use a
lower-degree factor `Q=V^(a-1)(V+W)`, giving `Phi((V+W)^2)=1`; when `a=1`,
the nonzero triple factor instead gives `Phi((V+W)^3)=4`.  The symmetric
argument handles a triple-factor boundary.  Together with irreducibility on
the monic chart, this shows that `X` is the projective closure of the
incidence hypersurface.  The only unstable exception is `(a,b)=(0,1)`: on
`P^1` with coordinates `[s:t]`,

\[
 \Phi((sV+tW)^3)=t^2(3s+t),                                  \tag{26}
\]

so the global divisor contains the base point with multiplicity two.  After
saturation the remaining point has degree one; equivalently, the unsaturated
moving-intersection calculation below subtracts the boundary length two.

The admissible monic open is dense in `X`, and the component-normalization
theorem makes (24) generically one-to-one there.  Thus its moving top
projective degree is the degree of `C_(a,b)`.

### Weighted local complete intersection

Work in the affine chart at `p`:

\[
 Q=1+\sum_{i=1}^a q_iW^i,
 \qquad
 R=1+\sum_{i=1}^b r_iW^i,
 \qquad
 M=Q^2R^3=\sum_jm_jW^j,
\]

and give `q_i,r_i` weight `i`.  Then `m_j` is weighted homogeneous of
degree `j`.  In the local polynomial ring

\[
 A=k[q_1,\ldots,q_a,r_1,\ldots,r_b]_{\mathfrak m}
\]

consider

\[
 J=(m_2,m_3,\ldots,m_{\ell+1}).                              \tag{27}
\]

This ideal is `m`-primary.  If its generators vanish, write

\[
 M=1+cW+O(W^{\ell+2})
\]

and put

\[
 B=(1+cW)(2Q'R+3QR')-cQR.
\]

The identity

\[
 QR^2B=(1+cW)M'-cM                                           \tag{28}
\]

shows that `B` vanishes to order at least `ell+1`, while `deg B<=ell`.
Hence `B=0` and `(1+cW)M'=cM`, which forces `M=1+cW`.  Unique factorization
then gives `c=0` and `Q=R=1`, because a nonconstant factor of `Q^2R^3`
cannot occur with multiplicity one.  Thus (27) has only the origin as a
geometric zero.

The ring `A` is Cohen--Macaulay of dimension `ell`, so the `ell` generators
in (27) form a regular sequence.  Its weighted Hilbert series gives

\[
 \begin{aligned}
 \operatorname{length}(A/J)
 &=\lim_{z\to1}
 {\prod_{j=2}^{\ell+1}(1-z^j)
  \over
  \prod_{i=1}^a(1-z^i)\prod_{i=1}^b(1-z^i)}\\
 &={2\cdot3\cdots(\ell+1)\over a!b!}
  ={(\ell+1)!\over a!b!}
  =(\ell+1)\binom{\ell}{a}.                                 \tag{29}
 \end{aligned}
\]

### Identification with the Segre correction

Let `I=(m_2,...,m_N)A` be the retained-coefficient base ideal.  On the local
hypersurface `A/(Phi)`, the zero-dimensional Segre contribution is the
Samuel multiplicity of `I A/(Phi)`.  Choose `ell-1` general linear
combinations `g_1,...,g_(ell-1)` of the retained coefficients.  Then

\[
 e\bigl(I A/(\Phi)\bigr)
 =\operatorname{length}
 A/(\Phi,g_1,\ldots,g_{\ell-1}).                              \tag{30}
\]

For a general choice, the coefficient matrix of these `ell` equations in
the columns `m_2,...,m_(ell+1)` is invertible.  Row reduction replaces them
by

\[
 h_j=m_j+\sum_{k\ge\ell+2}c_{jk}m_k,
 \qquad2\le j\le\ell+1.                                     \tag{31}
\]

Thus `in_w(h_j)=m_j`.  Since these initial forms are a regular sequence, the
filtered complete-intersection lemma identifies the associated graded
quotient with `A/J`.  Concretely, homogenize each `h_j` with a Rees parameter
for the weight filtration.  The Koszul complex of the special fiber is exact
because `m_2,...,m_(ell+1)` is a regular sequence; the same filtered Koszul
argument gives a flat finite family, so the special and local general fibers
have the same length.  Equations (29)--(31) therefore give the local Segre
correction

\[
 \boxed{e_p=(\ell+1)\binom{\ell}{a}.}                         \tag{32}
\]

Finally,

\[
 \int_Xc_1(L)^{\ell-1}
 =\int_S(2u+3v)^\ell
 =\binom{\ell}{a}2^a3^b.
\]

Subtracting (32) proves (22).  Summing the raw terms gives (19), while

\[
 \sum_{a,b\ge0}(\ell+1)\binom{\ell}{a}
 (tx^2)^a(tx^3)^b
 ={1\over(1-tx^2-tx^3)^2},
\]

which proves (23).  The checks `(a,b)=(2,0),(1,1),(3,0)` give degrees
`1,6,4`, respectively, agreeing with the dense quartic component, the
explicit quintic exceptional equation, and the explicit degree-six quartic
surface.

### Full graph multidegrees and Chow class

The same zero-dimensionality makes every graph multidegree immediate.  Let
`Gamma_(a,b)` be the closure of the graph of (24), let `h` be the target
hyperplane class, and define

\[
 \gamma_{p,q,r}^{a,b}
 =\int_{\Gamma_{a,b}}u^pv^qh^r,
 \qquad p+q+r=\ell-1.                                       \tag{33}
\]

If `p+q>0`, choose the source hyperplanes generally so that their intersection
does not contain `p`.  Their inverse image in the graph is then disjoint from
the exceptional locus over the base point.  On this intersection the graph
is the graph of a morphism, so the target hyperplanes pull back to `L` and
there is no Segre correction.  Only `(p,q,r)=(0,0,ell-1)` sees the base point.
Consequently

\[
 \boxed{
 \gamma_{p,q,r}^{a,b}
 =\binom{r+1}{a-p}2^{a-p}3^{b-q}
 -\mathbf1_{p=0}\mathbf1_{q=0}
  (\ell+1)\binom{\ell}{a}.}                                 \tag{34}
\]

As usual, the binomial coefficient is zero outside its natural range.  Thus
all the raw source multidegrees (20) survive unchanged except the top target
degree, where (34) becomes (22).

There is a compact generating series for the entire collection.  With
independent variables `A,B` marking `(a,b)`, `S,T` marking source
hyperplanes, and `Z` marking target hyperplanes, put

\[
 \mathscr G(A,B;S,T,Z)
 =\sum\gamma_{p,q,r}^{a,b}A^aB^bS^pT^qZ^r,
\]

where the sum is over `a,b,p,q,r>=0` with `p+q+r=a+b-1`.  Then

\[
 \boxed{
 \mathscr G
 ={2A+3B\over
   (1-AS)(1-BT)(1-Z(2A+3B))}
 -{1\over Z}\left({1\over(1-Z(A+B))^2}-1\right).}            \tag{35}
\]

The first term sums (20): write `a=p+alpha`, `b=q+beta`, so
`alpha+beta=r+1`.  The second term is the unique-base-point correction (32),
supported at `p=q=0`.

Finally, put `m=N-3`.  Since the Chow ring of
`S x P^m` is generated by `u,v,h`, these multidegrees determine the graph
class itself:

\[
 \boxed{
 [\Gamma_{a,b}]
 =\sum_{p+q+r=\ell-1}
 \gamma_{p,q,r}^{a,b}\,
 u^{a-p}v^{b-q}h^{m-r}
 \quad\text{in }A^{m+1}(S\times\mathbb P^m).}                 \tag{36}
\]

Thus the retained-coefficient projective degrees, all graph multidegrees,
and the graph Chow classes are complete.  Minimal pairwise intersection
weights are computed in the partition-lattice theorem; nonminimal collision
blocks and complete multi-sheet conductor equalizers still require their
normal cones and conductors.

## 7. Universality of the local correction

The correction above is not special to the cusp atoms two and three.  Let
`e_1,...,e_s>=2` be labelled multiplicity atoms, let `deg Q_nu=a_nu`, and put

\[
 M=\prod_{\nu=1}^sQ_\nu^{e_\nu},
 \qquad
 \ell=\sum_{\nu=1}^sa_\nu.
\]

In the chart at `Q_nu=1`, give the coefficient of `W^i` in every `Q_nu`
weight `i`, and again write `M=sum_j m_jW^j`.  Then

\[
 J_{\boldsymbol e}=(m_2,\ldots,m_{\ell+1})                  \tag{37}
\]

is primary to the coefficient maximal ideal and is a weighted complete
intersection.  To see this, set

\[
 P=\prod_\nu Q_\nu,
 \qquad
 C=\sum_\nu e_\nu Q_\nu'\prod_{\mu\ne\nu}Q_\mu.
\]

If (37) vanishes and `M=1+cW+O(W^(ell+2))`, define

\[
 B=(1+cW)C-cP.
\]

Since

\[
 \left(\prod_\nu Q_\nu^{e_\nu-1}\right)B
 =(1+cW)M'-cM,                                               \tag{38}
\]

the same order--degree comparison gives `B=0`, hence
`(1+cW)M'=cM` and `M=1+cW`.  A product of powers with all exponents at least
two cannot have a simple root, so `c=0` and every `Q_nu=1`.  Therefore

\[
 \boxed{
 \operatorname{length}(A/J_{\boldsymbol e})
 ={(\ell+1)!\over\prod_\nu a_\nu!}
 =(\ell+1)\binom{\ell}{a_1,\ldots,a_s}.}                     \tag{39}
\]

In particular, the local Segre correction depends only on the numbers of
support points of each labelled type, not on their multiplicities.  The
multiplicities enter only through the raw polarization
`O(e_1,...,e_s)`.

This gives a conditional global formula useful beyond the present 2/3
geometry.  If the corresponding projective incidence divisor has no
additional boundary component, its total moving top degree is

\[
 \binom{\ell}{a_1,\ldots,a_s}
 \left(\prod_\nu e_\nu^{a_\nu}-(\ell+1)\right).              \tag{40}
\]

If that divisor is irreducible and its retained-coefficient projection has
generic degree `delta`, then (40) equals `delta` times the degree of its
image.  Thus irreducibility and generic degree are the only global inputs not
provided by the local lemma.  The one-support-point compactification has the
same harmless exception as (26): `Phi((sV+tW)^e)` contains `t^2`, and
saturation or subtraction removes exactly the universal correction two.

For any fixed finite list of labelled atoms, the total moving-degree series
is consequently

\[
 \boxed{
 {1\over1-t\sum_\nu e_\nu x^{e_\nu}}
 -{1\over\left(1-t\sum_\nu x^{e_\nu}\right)^2}.}             \tag{41}
\]

Taking `(e_1,e_2)=(2,3)` recovers (22)--(23).

## Executable certificate

Run

```bash
.venv/bin/python scripts/verify_nonsurjective_enumerative_geometry.py
```

The script checks the quasipolynomial and all rational generating functions,
the component staircases, the interval-intersection theorem through degree
twenty-two, the binomial normalization fibers, the coincident-root Chow and
multidegree formulas, the weighted local lengths through support four, and
the independent degree-six sextic-curve regression.  It also checks the full
graph-multidegree generating series through support five and several local
systems with multiplicity atoms beyond two and three.
