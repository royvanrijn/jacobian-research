# Consecutive-factor contact reduction

This note isolates the remaining geometry of the normalized consecutive
factor spaces. It gives an exact motivic reduction for every hyperplane and
records finite-field evidence in the first unresolved degree. The reduction
is a classification algorithm, not yet a closed formula for all contact
moduli.

Work over a field `k`. Put

\[
 V_d=\operatorname{Sym}^d(k^2),\qquad n=2p+1,
\]

and let `ell in V_n^*` be nonzero. On
`P(V_p) x P(V_(p+1))` write

\[
 C_p=\{(A,B):\gcd(A,B)=1\},\qquad
 S_{p,\ell}=\{(A,B)\in C_p:\ell(AB)=0\},
\]

and `U_(p,ell)=C_p setminus S_(p,ell)`. For consecutive degrees, unique
resultant--hyperplane normalization identifies this open with the normalized
affine source `Res(A,B)=ell(AB)=1`.

## 1. The universal coprime class

Let `L=[A^1]`. For every `a,b>0`,

\[
 \boxed{[\operatorname{Cop}_{a,b}]
   =\mathbb L^{a+b-1}(\mathbb L+1).}                 \tag{1}
\]

Indeed, stratify `P^a x P^b` by the degree `d` of the projective gcd. The
`d`-th stratum is

\[
 \mathbb P^d\times\operatorname{Cop}_{a-d,b-d},
\]

with endpoint class `[P^e]` when the other residual degree is zero.
Induction and the geometric-series identity give (1). In particular,

\[
 [C_p]=\mathbb L^{2p}(\mathbb L+1).                  \tag{2}
\]

Consequently a necessary condition for an affine-space exception is

\[
 \boxed{[S_{p,\ell}]=\mathbb L^{2p}.}                \tag{3}
\]

This is the most economical motivic test: only the hyperplane section of the
coprime locus has to be computed.

## 2. Catalecticant term and exact gcd recursion

The hyperplane equation defines the rectangular catalecticant

\[
 \kappa_{p,\ell}:V_p\longrightarrow V_{p+1}^*,
 \qquad A\longmapsto(B\mapsto\ell(AB)).              \tag{4}
\]

Let `K_(p,ell)=P(ker(kappa_(p,ell)))`, with the empty projective space having
class zero. The complete bilinear hyperplane has class

\[
 \boxed{[E_{p,\ell}]=[\mathbb P^p]^2
       +\mathbb L^{p+1}[K_{p,\ell}].}                \tag{5}
\]

For `1<=d<=p`, define

\[
 Z_{d,\ell}=\{(D,A',B'):\ D\in\mathbb P^d,
  (A',B')\in\operatorname{Cop}_{p-d,p+1-d},
  \ell(D^2A'B')=0\}.                                \tag{6}
\]

Every pair on `E_(p,ell)` has a unique projective gcd. Its exact gcd-degree
strata give

\[
 \boxed{[S_{p,\ell}]=[\mathbb P^p]^2
 +\mathbb L^{p+1}[K_{p,\ell}]
 -\sum_{d=1}^p[Z_{d,\ell}].}                        \tag{7}
\]

Combining (2) and (7),

\[
 \boxed{[U_{p,\ell}]=\mathbb L^{2p}(\mathbb L+1)
 -[\mathbb P^p]^2-\mathbb L^{p+1}[K_{p,\ell}]
 +\sum_{d=1}^p[Z_{d,\ell}].}                       \tag{8}
\]

Equations (4)--(8) are valid in the Grothendieck ring and, with brackets
replaced by point counts, over every finite field. They retain all contact
moduli: on the fiber over `D`, (6) is the lower-degree hyperplane defined by

\[
 \ell_D(F)=\ell(D^2F).                               \tag{9}
\]

Thus the recursion does not collapse distinct cross-ratio strata.

## 3. Contact data and the classification algorithm

Restricting `ell` to the rational normal curve gives

\[
 h_\ell(L)=\ell(L^{2p+1}).                           \tag{10}
\]

Its multiplicity partition is the first invariant. Once its support has at
least four points, the configuration in `M_(0,r)/S_lambda` is also part of
the stratum. Formula (9) shows why this datum cannot be discarded: the
family of lower catalecticants obtained while `D` varies depends on the
actual configuration, not only on the multiplicities in (10).

The finite classification procedure is:

1. fix a multiplicity partition of (10) and its configuration-moduli stratum;
2. compute the rank strata of (4);
3. stratify `P^d` by the contact type of the contracted form (9);
4. evaluate (7), beginning with the classified `(1,2)` row;
5. retain only strata satisfying (3), and use point counts before attempting
   coordinates or locally nilpotent derivations.

## 4. The quadratic--cubic row

For `p=2`, (7) becomes

\[
 [S_{2,\ell}]=[\mathbb P^2]^2+\mathbb L^3[K_{2,\ell}]
                -[Z_{1,\ell}]-[Z_{2,\ell}],         \tag{11}
\]

where `Z_(1,ell)` is a `P^1`-family of linear--quadratic coprime hyperplane
sections and

\[
 Z_{2,\ell}=\{(D,L)\in\mathbb P^2\times\mathbb P^1:
                    \ell(D^2L)=0\}.                 \tag{12}
\]

Thus every quintic hyperplane reduces to rank strata of a `3 x 4`
catalecticant, the known cubic three-orbit calculation, and (12). This is
substantially smaller than direct elimination in seven affine coefficients.

For `ell(AB)=[AB]_(T^4S)`, (11) recovers

\[
 [S_{2,\ell}]=\mathbb L^4+\mathbb L^3,
 \qquad [U_{2,\ell}]=\mathbb L^5-\mathbb L^3,        \tag{13}
\]

so the natural tangent slice fails (3).

An exhaustive point-count sweep gives an additional experimental fact. For
every degree-five hyperplane over `F_3`, `F_5`, and `F_7` (respectively 364,
3906, and 19608 hyperplanes), no complement has `q^5` points. The number of
distinct complement counts is respectively 9, 11, and 13. In particular,
the count is not determined by the multiplicity partition alone; examples
with five simple rational contact points change with their cross-ratio
configuration.

This sweep is evidence, not a characteristic-zero nonexistence proof. A
proof that `(1,2)` is the unique affine-space exception still requires a
uniform evaluation or obstruction for the moduli-dependent terms in (11),
and then for (7) in arbitrary `p`.

## 5. Complete classification of the quadratic--cubic row

The recursion (11) admits a closed form. For `D in P^1`, put

\[
 c_D(L)=\ell(D^2L^3).
\]

Let `T` be the locus where `c_D` is singular or zero, `O_3` the locus where
it is a nonzero cube, and `O_0` the locus where it is zero. Finally put

\[
 Q=\{D\in\mathbb P(V_2):\ell(D^2V_1)=0\},\qquad
 K=\mathbb P(\ker\kappa_{2,\ell}).                   \tag{14}
\]

All loci have their reduced structure. The cubic `(1,2)` classification gives

\[
 [Z_{1,\ell}]=\mathbb L^2[\mathbb P^1]
 +\mathbb L([\mathbb P^1]-[T])
 +\mathbb L^2[O_3]+\mathbb L^3[O_0].                \tag{15}
\]

Indeed, the coprime cubic hyperplane section has class `L^2+L`, `L^2`,
`2L^2`, or `L^3+L^2` according as its type is `(1,1,1)`, `(2,1)`, `(3)`,
or the functional is zero. Projection of (12) to `P(V_2)` similarly gives

\[
 [Z_{2,\ell}]=[\mathbb P^2]+\mathbb L[Q].           \tag{16}
\]

Substitution in (11) yields

\[
 \boxed{[U_{2,\ell}]=\mathbb L^5-\mathbb L^3
 +\mathbb L([Q]-[T])+\mathbb L^2[O_3]
 +\mathbb L^3([O_0]-[K]).}                          \tag{17}
\]

### Theorem 5.1

Over a characteristic-zero field, no normalized quadratic--cubic
hyperplane complement is affine five-space:

\[
 \boxed{U_{2,\ell}\not\simeq\mathbb A^5
        \quad\text{for every }0\ne\ell\in V_5^*.}  \tag{18}
\]

### Proof

Apply the Hodge--Deligne realization to (17), after base change to `C`, and
write `t=uv`. The catalecticant (4) has rank one, two, or three.

If its rank is three, `K` is empty. Since `O_0` is the inverse image of `K`
under the Veronese map `D -> D^2`, it is empty as well. The loci `Q` and `T`
have dimension at most one. The locus `O_3` is finite: if the quadratic map

\[
 \mathbb P^1\dashrightarrow\mathbb P(V_3),\qquad
 D\longmapsto[c_D]
\]

contained a dense set in the twisted cubic, its image would be constant.
Indeed, a nonconstant map of degree `e` to the twisted cubic pulls back its
hyperplane bundle with degree `3e`, whereas the displayed coordinates have
degree at most two after removal of their common factor. Coefficient
comparison in the constant case makes the catalecticant rank one. Thus every
correction term in (17) has Hodge degree at most two, and the coefficient of
`t^3` remains `-1`.

If the catalecticant has rank two, `K` is one point and `O_0` is either empty
or its single Veronese preimage. The same degree argument makes `O_3` finite.
Hence the coefficient of `t^3` in (17) is `-2` or `-1`, again not zero.

If the catalecticant has rank one, its Hankel minors say that `ell` is a pure
fifth power in the dual variables. Projective `SL_2` equivalence reduces to
an osculating coefficient hyperplane. Here

\[
 K\simeq Q\simeq T\simeq\mathbb P^1,\qquad
 O_3\simeq\mathbb A^1,qquad O_0\simeq\mathrm{pt},
\]

and (17) gives

\[
 [U_{2,\ell}]=\mathbb L^5-\mathbb L^4.              \tag{19}
\]

In all three ranks the Hodge--Deligne polynomial differs from `t^5`, proving
(18). \(\square\)

Thus `(1,2)` is the only affine-space exception among the first two
consecutive rows. Formula (17) also explains the finite-field variation:
Frobenius acts on `Q,T,O_3,O_0`, so equal contact partitions need not have
equal counts.

## 6. Topology and LND tests

For the natural `p=2` tangent slice, the proposed next cohomology test is
already settled:

\[
 H^2(U_{2,\ell},\mathbb Q)=0,\qquad
 H^3(U_{2,\ell},\mathbb Q)=\mathbb Q(-2).
\]

These groups imitate `A^2 x SL_2`, not `A^5`; the motivic class (13) already
excludes affine space. Searching for two commuting locally nilpotent
derivations with slices is therefore useful only as a product or
automorphism-group question for this failed candidate. For discovering an
affine-space exception, the efficient order is (3), finite-field counts,
topology, and only then LND coordinates.
