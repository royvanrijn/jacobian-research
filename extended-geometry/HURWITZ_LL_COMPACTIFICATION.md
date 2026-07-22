# The Hurwitz--Lyashko--Looijenga compactification of weighted seeds

This note identifies the discriminant pencil with a universal branch-value
incidence, constructs a proper admissible-cover closure of the normalized
seed locus, and explains how rerooting, collided primitive roots, the
conductor, and the reconstruction open fit into that closure.

Work over a characteristic-zero field.  Let `H` be a normalized admissible
degree-`N` seed:

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,
\]

and put

\[
 P_{H,s}(W)=H(W)-sW.
\]

## 1. The discriminant is an LL spectral incidence

For a degree-`N` polynomial `P`, its Lyashko--Looijenga divisor is the
unordered critical-value divisor

\[
 \operatorname{LL}(P)
 =\sum_{dP(r)=0}(e_r-1)[P(r)]
 \in \operatorname{Sym}^{N-1}(\mathbb A^1).
\]

The usual LL morphism sends `P` to this divisor.  The repository curve is not
the LL morphism itself; it is the pullback of the universal point-in-branch-
divisor incidence along the affine line

\[
 \ell_H:\mathbb A^1_s\longrightarrow\operatorname{Poly}_N,
 \qquad s\longmapsto P_{H,s}.                          \tag{1.1}
\]

Indeed, a point `r` is critical for `P_(H,s)` precisely when

\[
 s=H'(r),
\]

and its critical value is

\[
 P_{H,s}(r)=H(r)-rH'(r)=-t.
\]

Consequently

\[
 D_H=\{(s,t):-t\in\operatorname{Supp}
                   \operatorname{LL}(P_{H,s})\},       \tag{1.2}
\]

and its critical-point incidence normalization is

\[
 r\longmapsto\bigl(H'(r),rH'(r)-H(r)\bigr).           \tag{1.3}
\]

This translates the repository strata into standard LL language:

- `H''(r)=0` is the caustic, where a critical point degenerates;
- `nu_H(r)=nu_H(u)` with `r!=u` is the Maxwell locus, where distinct
  critical points have the same critical value;
- a multiple primitive root is a ramification point of `P_(H,0)=H` over the
  marked value zero; and
- the conductor is the finite gluing algebra of the normalized LL incidence.

Degrees of LL maps on caustic, Maxwell, and more general Hurwitz strata are
the setting of Lando--Zvonkine,
[Counting ramified coverings and intersection theory on spaces of rational
functions I](https://arxiv.org/abs/math/0303218), and Zvonkine,
[Part II](https://doi.org/10.17323/1609-4514-2007-7-1-135-162).

## 2. The generic rerooting cover is a marked Hurwitz forgetful map

On the exact-double, distinct-extra-root locus one can write uniquely

\[
 H(W)=
 -\frac{W^2(W-1)\prod_{i=1}^{N-3}(W-\rho_i)}
        {\prod_{i=1}^{N-3}(1-\rho_i)}.                 \tag{2.1}
\]

Thus the normalized seed open is

\[
 \mathscr S_N^\circ\simeq
 \operatorname{Conf}_{N-3}
  (\mathbb A^1\setminus\{0,1\})/S_{N-3}.              \tag{2.2}
\]

Equivalently, it parametrizes polynomial Hurwitz covers

\[
 H:(\mathbb P^1;p_0,p_1,p_\infty)
     \longrightarrow(\mathbb P^1;q_0,q_\infty)
\]

with total ramification over `q_infinity`, a ramified point `p_0` over
`q_0`, and a marked unramified point `p_1` over `q_0`.  The three source
marks give the affine coordinate `W` with values `0,1,infinity`; the target
coordinate and its scale are fixed by `q_0=0`, `q_infinity=infinity`, and
`H'(p_1)=-1`.

Forgetting `p_1` gives the rerooting groupoid.  Generically the zero fiber
has profile

\[
 (2,1^{N-2}),
\]

so there are exactly `N-2` choices for `p_1`.  Choosing the point `a` gives

\[
 \kappa_a=-\frac1{aH'(a)},\qquad
 G_a(w)=\kappa_aH(aw),                                \tag{2.3}
\]

which is the repository rerooting identity.  The full reconstruction open
retains the marked regular component and therefore rigidifies this
forgetful map; the coarser decorated normalization forgets it and has generic
degree `N-2`.

## 3. Root-collision and admissible-cover compactifications

There are two compatible levels of compactification.

### 3.1 Collision-separating source compactification

Label the `N-3` extra roots temporarily.  Adding the three heavy marks
`0,1,infinity` identifies the labelled open with an open of `M_(0,N)`.
Define

\[
 \overline{\mathscr S}^{\mathrm{root}}_N
 =\overline M_{0,N}/S_{N-3}.                           \tag{3.1}
\]

This maximal model replaces every collision by a rational bubble and retains
the complete collision tree.  Hassett weights on the extra roots contract
chosen bubbles and allow collisions of bounded total weight on a smooth
component.  Varying the weights therefore interpolates between the fully
separated tree and the coefficient-space closure with repeated roots.  See
Hassett,
[Moduli spaces of weighted pointed stable curves](https://arxiv.org/abs/math/0205009).

### 3.2 Marked admissible-cover compactification

Let `\mathscr H_N^{\mathrm{mark}}` be the Hurwitz stack of the framed
polynomial covers
above, with all remaining ramification unordered.  Embed
`\mathscr S_N^\circ` in the stack of marked admissible covers with profile
`(N)`
over `q_infinity`, and define

\[
 \overline{\mathscr S}^{\mathrm{adm}}_N
 =\operatorname{Norm}
   \overline{\mathscr S_N^\circ}^{\,\mathrm{Adm}}.     \tag{3.2}
\]

This is a proper Deligne--Mumford stack: the admissible-cover stack is proper,
and normalization of the finite-type closure is finite.  Its points retain
the source and target bubble trees, equal ramification indices on the two
branches above every target node, the marked zero cluster, and the marked
regular root-one branch.  The normalization has the twisted-cover modular
interpretation of Abramovich--Corti--Vistoli,
[Twisted bundles and admissible covers](https://arxiv.org/abs/math/0106211).
Deopurkar's
[Compactifications of Hurwitz spaces](https://arxiv.org/abs/1206.4535)
provides intermediate models in which branch points may collide to a
prescribed extent; these are the Hurwitz analogue of the Hassett contractions
above.

The morphism from (3.1) to (3.2) is obtained by stable reduction of the
polynomial map determined by (2.1).  It remembers root-collision directions
on the source, while the admissible target also separates colliding critical
values.  Taking its normalization is essential: it retains the roots of the
smoothing parameters which appear as boundary inertia.

## 4. Local collision models and the compactified rerooting groupoid

Suppose `mu` simple zero-fiber points collide at a nonzero root `rho`.  A
transverse labelled model is

\[
 H_\epsilon(W)=U(W)\bigl((W-\rho)^\mu-\epsilon\bigr),
 \qquad U(\rho)\ne0.                                  \tag{4.1}
\]

After normalization write `epsilon=delta^mu`.  The `mu` choices of a root are
cycled by `delta -> zeta_mu delta`.  At the limit the marked rerooting factor
has

\[
 H_\epsilon'(a)\asymp\delta^{\mu-1},
 \qquad \kappa_a\asymp\delta^{1-\mu};                 \tag{4.2}
\]

so the rerooted seed leaves the affine normalized-seed chart.  The
admissible-cover bubble absorbs this divergent scale.  The completed
forgetful map has inertia `mu`, exactly matching the repository boundary
chart

\[
 C=\delta^\mu,\qquad W=\rho+\kappa\delta.              \tag{4.3}
\]

If extra roots collide with zero, the contracted coefficient-space seed has
`ord_0(H)>2`; the source-stable model retains the colliding points on a bubble.
For a zero of multiplicity `m`, the Newton polygon of

\[
 hW^m-BCW+AC^2                                      \tag{4.4}
\]

has two pieces.  For `m>=3`, one slope-one root remains in the reconstruction
open and `m-1` roots form the ramified boundary cluster; for `m=2`, both roots
belong to the finite quadratic chart.  This is precisely the affine/boundary
decomposition used in the generalized faithfulness theorem.

## 5. Extending the LL incidence and the conductor

Fantechi--Pandharipande construct a branch divisor in flat families and a
branch morphism from stable maps to a symmetric product:

\[
 \operatorname{br}:\overline{\mathcal M}
 \longrightarrow\operatorname{Sym}^{2N-2}(\mathbb P^1).
\]

For the polynomial profile, `(N-1)[infinity]` is the fixed branch
contribution.  Removing it leaves the finite LL divisor of degree `N-1` used
in (1.2).

See
[Stable maps and branch divisors](https://arxiv.org/abs/math/9905104).
Pulling the universal point-in-branch-divisor incidence back to
`\overline{\mathscr S}^{\mathrm{adm}}_N` extends (1.2).  Normalize that
incidence and
retain its conductor square

\[
 \begin{array}{ccc}
 \widetilde C&\longrightarrow&\widetilde D\\
 \downarrow&&\downarrow\\
 C&\longrightarrow&D.
 \end{array}                                           \tag{5.1}
\]

On the ordinary locus this square reduces to the cusp branches and paired
node branches already extracted in the repository.  At a collision over
`(s,t)=(0,0)`, it retains all ramification points over the common branch
value, not merely pairwise node data.  In curve language, (5.1) is the finite
normalization--conductor pushout; in admissible-cover language, it is the
algebraic shadow of the marked ramification points on the source bubble tree.

This gives the compactified invariant package

\[
 \boxed{
 (\text{marked admissible cover},\ \text{LL branch incidence},\
  \text{normalized conductor square},\ \text{reconstruction marking}) .}
                                                               \tag{5.2}
\]

Its open restriction is the repository full marked incidence cover.  Its
unmarked restriction is the finite rerooting groupoid, generically of degree
`N-2`; collision inertia and higher zero clusters are boundary phenomena of
the same proper stack rather than exceptional formulas.

## 6. The simultaneous-multicluster comparison theorem

The remaining comparison can be made without choosing a bubble chart.  For a
proper morphism `Y -> S`, write

\[
 \operatorname{St}^{\nu}(Y/S)
 =\operatorname{Norm}_S\operatorname{Spec}_S(f_*\mathcal O_Y)            \tag{6.1}
\]

for its normalized Stein factor.  The finiteness and relative-normalization
statement used here is the
[Noetherian Stein factorization theorem](https://stacks.math.columbia.edu/tag/03H0).

> **Theorem (formal admissible-cover comparison).**  Let `H` be any
> normalized admissible seed.  Pull the marked admissible-cover closure back
> along the pencil `s -> H-sW` and mark a point over the target value `-t`.
> Then:
>
> 1. the normalized Stein factor of the marked-point incidence is the
>    canonical finite root cover
>    \[
>      X_H=\operatorname{Spec}
>       k[s,t,W]/(H(W)-sW+t);
>    \]
> 2. the normalized Stein factor of its critical-point incidence is
>    `\widetilde D_H=Spec k[r]`, with finite map (1.3) onto `D_H`;
> 3. after contraction of the admissible source and target bubbles, every
>    completed boundary chart, the normalized LL branch incidence, and the
>    normalization--conductor square agree with the corresponding repository
>    constructions.
>
> These assertions remain true at simultaneous collisions of arbitrarily
> many distinct root clusters.

### 6.1 Uniqueness of the finite normal contraction

Away from the admissible-cover boundary, the marked-point incidence is
literally the equation `H(W)-sW+t=0`.  Its projection to the `(s,t)`-plane is
finite of degree `N`, and `X_H` is smooth because the derivative of its
defining equation with respect to `t` is one.  The normalized Stein factor
in (6.1) and `X_H` are therefore finite normal covers with the same generic
function algebra.  Both are the integral closure of `k[s,t]` in that algebra,
so they are canonically isomorphic.

The same argument on the critical marking gives the function field `k(r)`.
The finite normal model of the reduced scheme-theoretic image is uniquely
`Spec k[r]`, and its contraction is

\[
 r\longmapsto\bigl(H'(r),rH'(r)-H(r)\bigr).
\]

Thus taking an admissible-cover stable limit and then its normalized Stein
factor cannot manufacture a second formal model: normalization has already
made the contraction unique.

### 6.2 Completed root-cover charts

Let `alpha` be any primitive root of `H`, simple or multiple, put
`W=alpha+u`, and complete over `(s,t)=(0,0)`.  Eliminating `t` gives

\[
 \widehat{\mathcal O}_{X_H,(0,0,\alpha)}
 \simeq k[[s,u]],
 \qquad
 t=s(\alpha+u)-H(\alpha+u).                           \tag{6.2}
\]

Consequently the semilocal completion over the second boundary is the product
of the rings (6.2), one for each geometric root.  Hensel separation makes
distinct root clusters independent.  Formula (6.2) is also the contracted
local equation of the marked admissible cover, by Section 6.1.  The rational
reconstruction coordinates are the same functions in their common function
field.  Normality therefore gives the same pole divisors and the same maximal
regularity open on both finite models, including the separate affine and
polar valuations inside a higher root cluster.

### 6.3 The multicluster conductor

Suppose the distinct multiple roots are `rho_1,...,rho_q`, with
multiplicities `m_i>=2`.  Set

\[
 e_i=m_i-1,
 \qquad E=\sum_{i=1}^q e_i.
\]

Over a geometric residue field, write

\[
 H(\rho_i+u)=c_i u^{m_i}+O(u^{m_i+1}),\qquad c_i\ne0.
\]

On the corresponding branch of the LL incidence, put
`q_i=t-rho_i s`.  Directly from (1.3),

\[
 \operatorname{ord}_u(s)=e_i,
 \qquad
 \operatorname{ord}_u(q_i)=e_i+1.                    \tag{6.3}
\]

Thus the branch has semigroup `<e_i,e_i+1>` and intrinsic conductor exponent
`e_i(e_i-1)`.  Its tangent line is `t=rho_i s`.  Distinct root centers give
distinct tangent lines, so two different clusters have intersection number

\[
 I(D_i,D_j)=e_i e_j.                                  \tag{6.4}
\]

One can see the conductor addition directly from the repository adjunction
formula (3.2a): in the completed plane ring write the reduced equation as
`f=\prod_i f_i`.  On branch `i`, the derivative `f_t` is
`(f_i)_t\prod_{j\ne i}f_j`; the extra factors contribute precisely the
intersection numbers (6.4).  Therefore the conductor pulled back to the
completed normalization is

\[
 \boxed{
 \overline{\mathfrak c}_{(0,0)}
 =\bigoplus_{i=1}^q
   u_i^{\,e_i(E-1)}k[[u_i]].}                         \tag{6.5}
\]

This includes a higher zero cluster, any number of nonzero multiple-root
clusters, and every mixture of the two.  The completed local ring downstairs
is the image of

\[
 k[[s,t]]\longrightarrow\prod_{i=1}^q k[[u_i]]
\]

under the parametrizations (1.3), and its conductor square is obtained by
quotienting this map by (6.5).  But Section 6.1 identifies this very map with
the contracted admissible-cover LL incidence.  Hence the two conductor
squares coincide, not merely their reduced supports or numerical lengths.

Over a nonclosed field, the geometric factors and summands in (6.2)--(6.5)
are permuted by Galois.  Taking invariants groups them by the same residue
fields retained by the repository boundary extractor.  This proves the
comparison before and after geometric splitting.

The symbolic regression
`scripts/verify_multicluster_ll_comparison.py` checks (6.3)--(6.5) on
normalized seeds with two and three simultaneous clusters, including the
adjunction conductor obtained from the global implicit discriminant.  Thus
the local-to-global comparison left open after (5.2) is complete.
