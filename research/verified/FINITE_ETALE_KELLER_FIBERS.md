# Finite étale algebras as Keller fibers

## 1. Definition

Let \(K\) be a field. A nonzero finite étale \(K\)-algebra \(A\) is a
**Keller fiber** if there are a polynomial Keller map

\[
F:\mathbb A^m_K\longrightarrow\mathbb A^m_K
\]

and a point \(y\in\mathbb A^m(K)\) such that

\[
F^{-1}(y)\simeq\operatorname{Spec}A,
\qquad
\dim_KA=\operatorname{gdeg}(F).
\]

The equality says that the fiber is **full**: every generic inverse sheet is
present. A shorter special fiber can occur when sheets escape through the
nonproperness boundary; fullness excludes precisely that loss.

## 2. Effective Jacobian-one realization

Assume \(\operatorname{char}K\ne2\).  Let \(P\in K[T]\) be separable of
degree \(N\ge3\).  Write \(P^{[j]}\) for its \(j\)-th Hasse derivative,
defined by
\[
P(T+S)=\sum_{j\ge0}P^{[j]}(T)S^j,
\]
and choose \(a\in K\) with

\[
P^{[1]}(a)P^{[3]}(a)\ne0.
\]

and put

\[
G(S)=P(a+S)-P(a)=g_1S+\cdots+g_NS^N.
\]

Then \(g_1=P^{[1]}(a)\), \(g_3=P^{[3]}(a)\), and \(g_N\ne0\), so the
root-engineered quadratic gauge applies. Its raw map

\[
F_G=(\Pi,B,C):\mathbb A^3_K\longrightarrow\mathbb A^3_K
\]

has determinant \(-2\). Compose with the target automorphism

\[
L(\Pi,B,C)=\left(\Pi,-\frac{B}{2},C\right)
\]

and write \(\widetilde F_G=L\circ F_G\). Then

\[
\boxed{\det D\widetilde F_G=1},
\qquad
\operatorname{gdeg}(\widetilde F_G)=N.
\]

The distinguished target has second coordinate zero, so it is fixed by this
normalization:

\[
y_{P,a}=\left(1,0,-\frac{2P(a)}{P^{[1]}(a)}\right).
\]

At this target the full fiber is

\[
\boxed{
\widetilde F_G^{-1}(y_{P,a})\simeq
\operatorname{Spec}K[T]/(P).
}
\]

The construction is effective. With total degree,

\[
\deg\Pi\le7,
\qquad
\deg B\le6N+2,
\qquad
\deg C\le6N,
\]

so

\[
\boxed{\max_i\deg(\widetilde F_G)_i\le6N+2.}
\]

The estimates use \(\deg t=2\), \(\deg q\le5\), and

\[
\deg\bigl(t^2x^{k-2}q^k\bigr)\le6k+2,
\qquad
\deg\bigl((xq)^k\bigr)\le6k.
\]

## 3. Complete rank classification

Now assume that \(K\) has characteristic zero.
Every finite étale algebra over an infinite field is monogenic. After
extension to a separable closure, the algebra is a product of \(N\) copies of
the field; the primitive elements form the nonempty open set on which the
Vandermonde discriminant

\[
\prod_{i<j}(X_i-X_j)^2
\]

is nonzero. Thus every finite étale \(K\)-algebra of rank \(N\) is
\(K[T]/(P)\) for a squarefree degree-\(N\) polynomial.

Consequently every rank at least three occurs by the construction above.
Rank one is realized by the identity. Rank two is impossible by the
standalone degree-two lemma in the paper. Precisely, descend the coefficients
of a hypothetical map to a finitely generated field \(K_0/\mathbb Q\).
On a nonempty target open the map is finite locally free; its rank, and hence
its generic degree, is unchanged by every scalar extension. After choosing
\(K_0\hookrightarrow\mathbb C\), the complex function-field extension is
therefore still quadratic. It is separable and normal, so Campbell's
unnumbered theorem on p. 244 applies and makes the complexified map a
polynomial automorphism. The coordinate-ring homomorphism is then an
isomorphism already over \(K_0\) by faithful flatness, and hence over \(K\).
This contradicts degree two because an automorphism has function-field degree
one. Razar and Wright are later algebraic sources for the same Galois case.

This is an exclusion of the generic degree before a target point is chosen.
An empty fiber is not a nonzero Keller fiber under the definition above, and
a shorter special fiber of a map of another generic degree is not full.
Neither affects the classification. Hence the possible ranks of nonzero full
Keller fibers are exactly

\[
\boxed{1,3,4,5,\ldots}.
\]

## 4. Scheme-theoretic reconstruction

Return to an arbitrary field \(K\) of characteristic different from \(2\).
For a quadratic-gauge target \((\pi,b,c)\) with \(\pi\ne0\), set

\[
E(S)=G_\pi(S)-\frac{g_1}{2}(bS^2+c),
\qquad
R=K[S]/(E),
\qquad
s=S\bmod E.
\]

If \(E\) is separable, equivalently \(\gcd(E,E')=1\), Bézout gives
\(U,V\in K[S]\) with

\[
UE+VE'=1.
\]

Therefore \(E'(s)\) is a unit in \(R\), with explicit inverse \(V(s)\). Put

\[
d=\frac{E'(s)}{g_1},
\qquad
Q=b-\beta(\pi,s),
\qquad
\beta(\pi,S)=\frac{G_\pi'(S)/g_1-1-\pi S^2}{S},
\]

and reconstruct in \(R\)

\[
t=d^{-1},
\qquad
x=sd^{-1},
\qquad
y=Q-\pi s,
\qquad
q=\pi d,
\]

\[
z=d^2\left(q-\frac{g_1}{g_3}y^2(1+3t)\right).
\]

These are elements of the quotient ring, not pointwise rational functions.
Conversely, on the entire fiber ring,

\[
tq=\pi\in K^*,
\]

so \(t\) and \(q\) are units globally and

\[
S=x/t,
\qquad
Q=y+xq
\]

are global fiber-ring elements satisfying \(E(S)=0\). The two maps are inverse
on coordinate rings. Thus

\[
F_G^{-1}(\pi,b,c)\simeq\operatorname{Spec}K[S]/(E)
\]

scheme-theoretically.

## 5. Collision fiber and off-diagonal algebra

This section applies the collision-ring and marked cubic-sheet interface
recorded in the
[external collision-ideals audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md) to the
finite fibers constructed here.

Let \(A=K[T]/(P)\) be any of the squarefree rank-\(N\) fibers above. Base
change of the affine self-fiber product to its target point gives the ordered
collision fiber

\[
\operatorname{Spec}(A)\times_{\operatorname{Spec}K}\operatorname{Spec}(A)
=\operatorname{Spec}(A\otimes_KA).
\tag{5.1}
\]

Multiplication is the diagonal restriction

\[
\mu:A\otimes_KA\longrightarrow A,\qquad a\otimes b\longmapsto ab,
\tag{5.2}
\]

and the fiberwise collision obstruction is

\[
\operatorname{Obs}(A)=\ker\mu.
\tag{5.3}
\]

This is the base change of the global collision scheme to the distinguished
target. It records every ordered collision inside that complete fiber, but
does not by itself assert that a factor extends to a global affine component
of the map's collision scheme.

For the polynomial presentation, introduce two root coordinates \(r,u\) and
put

\[
D_P(r,u)=\frac{P(u)-P(r)}{u-r}.
\tag{5.4}
\]

Then

\[
A\otimes_KA
\simeq
\frac{K[r,u]}{(P(r),P(u))}.
\tag{5.5}
\]

Because \(P\) is squarefree, \(P'(r)\) is a unit modulo \(P(r)\). If
\(h(r)P'(r)\equiv1\pmod {P(r)}\), the element

\[
e_\Delta=h(r)D_P(r,u)
\tag{5.6}
\]

is an idempotent in (5.5). It is one on the diagonal and zero on the
off-diagonal factor. Thus the two ideals

\[
(P(r),u-r),\qquad (P(r),D_P(r,u))
\]

are comaximal and give the exact Chinese-remainder decomposition

\[
\boxed{
A\otimes_KA
\simeq
A\times A^{\mathrm{off}},\qquad
A^{\mathrm{off}}=
\frac{K[r,u]}{(P(r),D_P(r,u))}.
}
\tag{5.7}
\]

Writing \(q=1-e_\Delta\), multiplication in (5.2) is the first projection and

\[
\boxed{
\operatorname{Obs}(A)=(A\otimes_KA)q
\simeq A^{\mathrm{off}}.
}
\tag{5.8}
\]

The ranks are therefore

\[
\dim_K(A\otimes_KA)=N^2,\qquad
\dim_KA=N,\qquad
\dim_K\operatorname{Obs}(A)=N(N-1).
\tag{5.9}
\]

For every commutative \(K\)-algebra \(R\), the universal property gives

\[
\operatorname{Hom}_{K\text{-alg}}(A\otimes_KA,R)
\simeq
\operatorname{Hom}_{K\text{-alg}}(A,R)^2.
\tag{5.10}
\]

Thus the tensor-square local point count used in the Hasse argument is
literally the ordered collision-point count.

### Higher ordered-configuration fibers

The pair construction extends without choosing a polynomial presentation.
If \(E=\operatorname{Spec}(A)\), define

\[
\operatorname{Conf}_m(E/K)
=E^m\setminus\bigcup_{i<j}\Delta_{ij},
\qquad 1\le m\le N.
\tag{5.11}
\]

Because \(E/K\) is finite etale, every partial diagonal is clopen.  Hence
\(\operatorname{Conf}_m(E/K)\) is again finite etale, and after a separable
closure it is the set of injections from an ordered \(m\)-element set into
the \(N\) geometric points of \(E\).  Therefore

\[
\boxed{
\operatorname{rk}_K\operatorname{Conf}_m(E/K)
=N(N-1)\cdots(N-m+1)=\frac{N!}{(N-m)!}.
}
\tag{5.12}
\]

Let \(H\le S_N\) be the finite geometric monodromy image on those \(N\)
points.  The connected factors of \(\operatorname{Conf}_m(E/K)\) are
exactly the \(H\)-orbits on ordered injective \(m\)-tuples, with residual
fields fixed by the corresponding stabilizers.  In particular, at
\(m=N-1\) and \(m=N\) every stabilizer is trivial: every connected factor is
the full splitting field, although a special fiber can contain several
copies.  When \(H=S_N\), there is one orbit for every \(m\), its stabilizer
is \(S_{N-m}\), and the last two covers are single copies of the
\(S_N\)-normal closure.  The
[universal relative family](UNIVERSAL_RELATIVE_KELLER_MAP.md) realizes this
full-monodromy tower canonically.

### Cubic \(S_3\) fibers

If \(A=L\) is a separable nonnormal cubic field with normal closure \(N_L\),
then the off-diagonal algebra is a field and

\[
\boxed{
L\otimes_KL\simeq L\times N_L,\qquad
\operatorname{Gal}(N_L/K)\simeq S_3.
}
\tag{5.13}
\]

The first factor marks one root. The off-diagonal factor marks an ordered
second distinct root, and those two roots determine the third. Equation
(5.13) is a generic field statement when applied to the function field of a
map, but here it is an exact statement about the finite fiber algebra.

The same observation has a presentation-descent consequence.  For every
rank-three finite-etale algebra, including split products,
`Conf_2=Off_2` is canonically the full `S_3` frame torsor.  The
[rank-three collision-framed audit](RANK_THREE_COLLISION_DESCENT.md) uses
that torsor to construct the exact projective and quadratic-Tschirnhaus
cocycles between primitive presentations.  It lifts the projective
transition to the foundational factorization map after target localization;
global polynomial extension across its explicit denominator remains open.

In rank four, `Off_2` retains an `S_2` stabilizer and `Conf_3` is the full
`S_4` frame torsor.  The
[rank-four frame audit](RANK_FOUR_COLLISION_CROSS_RATIO.md) proves that a
primitive change
`u=q_0+q_1r+q_2r^2+q_3r^3` is projective on the four framed roots exactly
when

\[
 q_2^2-q_1q_3+q_2q_3e_1+q_3^2e_2=0.                 \tag{5.13a}
\]

Thus collision framing resolves the finite labeling ambiguity in both
ranks, but in rank four a genuine cross-ratio hypersurface remains before
any Keller-lift boundary is considered.

The
[all-rank projective-descent theorem](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md)
completes this calculation for every `N>=3`.  For primitive coordinates
`r,u in A`, put `W_r=<1,r>`.  Then

\[
 u=\frac{ar+b}{cr+d}\text{ for some }g\in PGL_2
 \quad\Longleftrightarrow\quad
 W_r\cap uW_r\ne0
 \quad\Longleftrightarrow\quad
 \operatorname{rank}(1,r,u,ru)\le3.                  \tag{5.13b}
\]

On `Conf_(N-1)`, which is the full `S_N` frame torsor, the last condition is
the vanishing of the `N-3` residuals left after projectively matching the
first three roots.  They cut out a smooth codimension-`N-3` locus on the
primitive-presentation open.  In normalized quotient coordinates they are
the `4`-by-`4` minors of one explicit `N`-by-`4` polynomial matrix.  Thus the
finite-etale collision tower removes all root-label ambiguity, while the
moduli of the embedded root configuration remain as genuine projective
obstructions.

The
[rank-four nonprojective Keller-lift theorem](RANK_FOUR_NONPROJECTIVE_KELLER_LIFT.md)
applies the first transverse direction to this exact fiber compiler.  It
shows that the quartic ground-field orbit has Kummer class

\[
 \frac{U'}U\in K^\times/K^{\times5},
\]

constructs a rational primitive nonprojective comparison with trivial
class, and reduces it to two regular fibers of the same map
`F_(-124416)`.  Those endpoints lie on the explicit finite-etale target
line

\[
 -\frac{(S-12)(S+12)(S^2+24S+108\lambda)}{3456},
\]

but the line's two constant sheets do not realize the required quadratic
collision-frame labels.  The desired framed motion nevertheless has an
exact divergence-free polynomial first-order lift and a unique lift to
every finite order.  At the line parameter `-4`, the inverse polynomial has
only two simple roots; the resulting four-to-two fiber-cardinality drop
proves that the straight target translation has no polynomial source lift.
The prime discriminant component has ordinary degree thirteen, so every
target self-equivalence through degree twelve is exactly in `mu_5`, and its
endpoint orbit also fails.  The exact logarithmic-boundary system and a
four-parameter Singular unit-ideal certificate further exclude endpoint
degrees thirteen through eighteen.  Global endpoint transport is now open
only from target degree nineteen onward, with the prescribed sheet
permutation.

For a product \(A=\prod_iL_i\), tensor distributivity refines the collision
algebra into the ordered blocks \(L_i\otimes_KL_j\). The cross blocks
\(i\ne j\) lie entirely in the obstruction, while each diagonal block
\(L_i\otimes_KL_i\) splits into its diagonal and internal off-diagonal
parts. This is the component surplus used by the second-moment argument.

For the explicit Hasse fiber

\[
A_5=\mathbb Q[T]/((T^3-19)(T^2+T+1))=L_3\times L_2,
\]

let \(N_6\) be the \(S_3\) normal closure of
\(\mathbb Q[T]/(T^3-19)\). Since
\(L_2=\mathbb Q(\sqrt{-3})\) is the quadratic resolvent of that cubic,

\[
L_3\otimes L_3\simeq L_3\times N_6,\qquad
L_3\otimes L_2\simeq N_6,\qquad
L_2\otimes L_2\simeq L_2\times L_2.
\]

Consequently

\[
\boxed{
A_5\otimes_{\mathbb Q}A_5
\simeq A_5\times
\bigl(N_6\times N_6\times N_6\times L_2\bigr),
}
\tag{5.14}
\]

with diagonal rank \(5\) and obstruction rank \(18+2=20\). This identifies
the precise collision algebra behind the minimal Hasse-failure calculation,
not only its first and second moments.

## 6. Scalar-extension compatibility

For every field extension \(K\hookrightarrow K'\) with
\(\operatorname{char}K\ne2\), the supplied \((P,a)\)-construction commutes
with coefficientwise scalar extension:

\[
\widetilde F_{G,K'}=\widetilde F_G\otimes_KK'.
\]

The distinguished fiber becomes

\[
(K[T]/(P))\otimes_KK'\simeq K'[T]/(P).
\]

Hence connectedness, field decomposition, real signatures, splitting fields,
Galois actions, and all local algebras are transported functorially.

In characteristic zero, automatic existence of an admissible translation
implies that every finite separable field extension \(L/K\) of degree at
least three occurs as a connected full fiber of a Jacobian-one map of
\(\mathbb A^3_K\), with coordinate degree at most \(6[L:K]+2\).

## 7. Explicit optimal Hasse fiber

Let

\[
P_5(T)=(T^3-19)(T^2+T+1)
\]

and

\[
G(T)=P_5(T)-P_5(0)=T^5+T^4+T^3-19T^2-19T.
\]

The normalized quadratic gauge has coefficient vector

\[
(g_1,g_2,g_3,g_4,g_5)=(-19,-19,1,1,1)
\]

and determinant \(-2\). If its target coordinates are \((\Pi,B_0,C_0)\), the
denominator-free displayed map is exactly

\[
(\Pi,19B_0,19C_0).
\]

Therefore its determinant is

\[
19^2(-2)=-722,
\]

and normalized target \((1,0,-2)\) becomes displayed target \((1,0,-38)\).
At the normalized target,

\[
G(T)-\frac{-19}{2}(-2)=G(T)-19=P_5(T).
\]

Hence

\[
F^{-1}(1,0,-38)\simeq
\operatorname{Spec}\mathbb Q[T]/((T^3-19)(T^2+T+1)).
\]

This fiber has points over every completion of \(\mathbb Q\) but no rational
point. A Dedekind-zeta first/second-moment argument excludes such finite
étale schemes in total degree at most four, so the minimum Hasse-failing
Keller-fiber rank is exactly five.  Equivalently, the first prime moment
counts global components, while the tensor square supplies the second
moment; a nontrivial field factor forces a strictly positive component
surplus and hence a contradiction.

## 8. Consequences for the earlier chain

1. **Absolute occurrence in characteristic zero is settled.** Every finite
   étale algebra of rank at least three occurs directly, with determinant one
   and an effective degree bound.  In characteristic different from two the
   same conclusion holds from any supplied separable presentation with an
   admissible translation.
2. Weighted tangent-admissibility remains useful for occurrence inside the
   specific weighted linear-pencil family, but it is not an absolute
   existence condition.
3. Quadratic-stabilized intersective transfer remains valid as a
   weighted-presentation theorem, but its degree-two overhead is unnecessary
   for general Keller-fiber existence.
4. The remaining arithmetic questions concern one fixed map, minimal
   coordinate complexity, stable-equivalence multiplicity, and additional
   restrictions on the ambient map.
5. The tensor-square second moment is the point count of the ordered
   collision fiber. Its strict component surplus is therefore a geometric
   off-diagonal-sheet surplus, not merely an analytic inequality.

The geometric-degree spectrum theorem now also has a second direct existence
proof: apply this realization to any squarefree polynomial of degree
\(N\ge3\).

## 9. Verification and formalization

Run

```bash
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

The exact checker verifies translated examples in degrees three, four, and
five, determinant-one normalization, the \(6N+2\) degree bound, Bézout
inversion, both quotient-ring reconstruction compositions, the explicit
quintic scaling, determinant \(-722\), target \(-38\), and the fixed-map
infinite Hasse-family identity. It also verifies the exact collision CRT for
the degree-three, degree-four, and degree-five examples; the ranks
\(N^2,N,N(N-1)\); three cubic \(S_3\) off-diagonal normal closures; and the
explicit decomposition
\(A_5\otimes A_5\simeq A_5\times(N_6^3\times L_2)\).
The full-\(S_N\) ordered-configuration ranks and stabilizers in (5.12) are
independently enumerated for every \(1\le m\le N\le8\) by
`scripts/verify_universal_relative_keller_map.py`.  The general
orbit--component statement is a written finite-etale/Galois-set argument,
not a conclusion of that bounded enumeration.

The two low-rank presentation audits are replayed by

```bash
.venv/bin/python scripts/verify_rank_three_collision_descent.py
.venv/bin/python scripts/verify_rank_four_collision_cross_ratio.py
.venv/bin/python scripts/verify_all_rank_collision_projective_descent.py
```

The first checks the full cubic frame, Tschirnhaus cocycle, target-localized
factorization transport, and global torus endpoint.  The second checks the
quartic frame stabilizers, fourth-root/cross-ratio defect, primitive
boundary, universal-target equation, and exact projective/nonprojective
witness cards.  The third checks the all-rank frame completion, coefficient
matrix, low-rank specializations, uniform witnesses, and `N-3` independent
projective residuals.  None of these checkers classifies arbitrary nonlinear
polynomial Keller equivalences.

A staged Lean project is stored in
[`formal/finite-etale-keller`](../formal/finite-etale-keller). It now
formalizes, in characteristic zero, the automatic three-variable realization,
full scheme-level fiber
reconstruction, naturality, finite étaleness, rank and geometric degree; the
explicit quintic Hasse certificate at every completion; and the entire
algebraic degree-four moment barrier. The module
[`CollisionFiber.lean`](../formal/finite-etale-keller/FiniteEtaleKeller/CollisionFiber.lean)
formalizes the self-fiber-product algebra, diagonal multiplication,
obstruction ideal, ordered-pair universal property, nonvanishing in rank
greater than one, and obstruction rank \(N^2-N\). For rank-minimality, only the
Dedekind-zeta Euler-product extraction of the first prime moment remains
outside Lean. The historical degree-two theorem remains a separate external
input until formalized.  The wider characteristic-not-two supplied theorem is
not yet formalized end to end.

The dated qualified novelty search is recorded in
[`papers/common-arithmetic-fibers/LITERATURE_AUDIT.md`](../papers/common-arithmetic-fibers/LITERATURE_AUDIT.md).
