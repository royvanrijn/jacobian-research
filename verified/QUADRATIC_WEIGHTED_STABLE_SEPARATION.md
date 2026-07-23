# Stable separation of the quadratic and weighted families

This note proves that the root-engineered quadratic-gauge family is not a
reparametrization of the weighted tangent family.  The comparison is made in
the canonical Zariski--Main normalization package, so it survives arbitrary
polynomial source and target automorphisms and arbitrary identity
stabilization.

Work over an algebraically closed field `k` of characteristic zero.  Fix
`N>=4`.

For the quadratic family take

\[
 G(S)=g_1S+\cdots+g_NS^N,\qquad g_1g_3g_N\ne0,
\]

and let `F_G^quad` be the map in
[the root-engineered quadratic-gauge theorem](../cancellation/ROOT_ENGINEERED_QUADRATIC_GAUGE.md).
Its inverse equation is

\[
 E(P,B,C;S)
 =G_P(S)-\frac{g_1}{2}(BS^2+C),                         \tag{1}
\]

where

\[
 G_P(S)=g_1S+P(g_2S^2+g_3S^3)
       +\sum_{j=4}^Ng_jP^jS^j.                          \tag{2}
\]

For the weighted family take a boundary-clean degree-`N` seed `H`: zero is
an exact double root and every nonzero root is simple.  Let `F_H^wt` be its
weighted Keller map, with the maximal regular-reconstruction open.  These
are precisely the hypotheses of the
[weighted intrinsic-boundary theorem](../papers/decorated-discriminant-normalization/main.tex).

## Stable separation theorem

No `F_G^quad` is stably polynomially left--right equivalent to a
boundary-clean `F_H^wt` of the same degree.

The same conclusion holds over every characteristic-zero ground field:
a stable equivalence over that field would remain one after extension to
its algebraic closure.

More precisely, both canonical normalization packages intrinsically select
an ordered pair of target boundary images `(Z_Delta,Z_0)`.  After deleting
`Z_0`, the normalization of the selected ramified target stratum is

\[
 \begin{array}{c|c|c}
 \text{family}&
 \operatorname{Norm}(Z_\Delta\setminus Z_0)&
 \operatorname{rank}\bigl(\mathcal O^\times/k^\times\bigr)\\ \hline
 \text{weighted}&\mathbb A^1\times\mathbb G_m&1\\
 \text{quadratic gauge}&\mathbb G_m^2&2.
 \end{array}                                             \tag{3}
\]

Polynomial stabilization adds an affine-space factor and changes neither
unit rank.  Hence the two intrinsic strata cannot become isomorphic.

The common `G_m` in (3) is important.  The plane pencils have normalization
`A^1` and `G_m`, respectively, but a fixed-coordinate plane slice is not
itself invariant under arbitrary target automorphisms.  The canonical
three-dimensional comparison is therefore

\[
 \mathbb A^1\times\mathbb G_m
 \quad\text{versus}\quad
 \mathbb G_m^2,                                         \tag{4}
\]

not merely `A^1` versus `G_m`.

## 1. The quadratic canonical boundary

Let

\[
 \overline X_G=\operatorname{Norm}_{\mathbb A^3_{P,B,C}}
 k(\mathbb A^3_{x,y,z})
\]

and let `partial_G` be the reduced complement of the affine source in
`\overline X_G`.  We first enumerate every codimension-one component of
`partial_G`, its target image, and its `(e,f)` labels.

### The ramified discriminant vertex

Assume `P!=0`.  A repeated inverse root `r` of (1) satisfies

\[
 B(r)=\frac{G_P'(r)}{g_1r},\qquad
 C(r)=\frac{2G_P(r)-rG_P'(r)}{g_1}.                     \tag{5}
\]

The critical root cannot be zero because `G_P'(0)=g_1`.  The identity

\[
 dC+r^2\,dB=0                                           \tag{6}
\]

holds at fixed `P`.  Put

\[
 K_0=k(P),\qquad
 K_\Delta=K_0(B(r),C(r))\subset K_0(r).
\]

Equation (6) gives `r^2 in K_Delta`.  Since `B(r)` has a pole of odd order
one at `r=0`, it does not lie in `K_0(r^2)`.  As
`[K_0(r):K_0(r^2)]=2`, this proves

\[
 K_\Delta=K_0(r).                                       \tag{7}
\]

Thus (5) is birational onto an irreducible reduced target divisor
`Z_Delta`.  At its generic point the root is exactly double.  With
`u=partial_S E=g_1D` as an upstairs uniformizer, the reduced local
discriminant equation satisfies

\[
 \delta_\Delta=\text{unit}\cdot u^2.                    \tag{8}
\]

There is consequently one boundary prime over `Z_Delta` with
`(e,f)=(2,1)`.  It is outside the affine source because `D=0` makes
`t=D^{-1}` singular.  The other `N-2` roots are simple and reconstruct
regular affine source points.  The local degree sum is

\[
 2+(N-2)=N.                                             \tag{9}
\]

### The second vertex `P=0`

Put

\[
 d=N-3,\qquad h=\gcd(d,2),\qquad e_0=\frac d h.          \tag{10}
\]

Over the generic DVR of `P=0`, with residue field `k(B,C)`, the coefficient
valuations of (1) have lower Newton polygon

\[
 (0,0)\longrightarrow(2,0)\longrightarrow(3,1)
 \longrightarrow(N,N).                                 \tag{11}
\]

All points `(j,j)` for `4<=j<N` lie strictly above the last segment.  The
three Newton blocks give the complete table

\[
\begin{array}{c|c|c|c}
 \text{root valuation}&\text{multiplicity}&(e,f)&\text{position}\\ \hline
 v(S)=0&2&(1,2)&q=0\ \text{affine}\\
 v(S)=-1&1&(1,1)&t=0\ \text{affine}\\
 v(S)=-(d+2)/d&d&
 h\text{ primes }(e_0,1)&\text{boundary}.
\end{array}                                             \tag{12}
\]

Here the first residual polynomial is

\[
 -\frac{g_1}{2}(BS^2-2S+C),                             \tag{13}
\]

which is irreducible over `k(B,C)` and is the generic degree-two map of the
affine divisor `q=0`.  On `q=0` one has

\[
 B=y,\qquad C=2S-BS^2,
\]

so (13) is exactly its residual inverse equation.

For the second block, write `S=P^{-1}Z`.  Its nonzero residual equation is

\[
 g_3Z-\frac{g_1}{2}B=0.                                 \tag{14}
\]

It is the affine divisor `t=0`: there `B=-2y`,
`q=(g_1/g_3)y^2`, and

\[
 PS=xq=\frac{g_1}{2g_3}B.
\]

Thus this branch also has `(e,f)=(1,1)`.

These are the only affine divisorial branches over `P=0`.  Indeed
`P=tq` in the factorial source ring, `t=1+xy` is irreducible, and `q` is
primitive and linear in `z`: its coefficients `t^2` and
`(g_1/g_3)y^2(1+3t)` are coprime.  Hence `q` is irreducible as well.

For the final block the endpoint equation is, up to a monomial factor,

\[
 g_3+g_NZ^d=0.                                         \tag{15}
\]

Equivalently, after inversion of the root coordinate, its binomial model is

\[
 W^d=\text{unit}\cdot P^{d+2}.
\]

Since `h=gcd(d,d+2)`, this has `h` geometric primes, each totally ramified
with index `e_0=d/h` and residue degree one.  The intermediate terms in (2)
lie strictly above the Newton segment, while (15) is separable, so Newton
lifting introduces no additional prime.

These last primes are genuinely outside the affine source.  Normalize one
of their valuations by

\[
 v(P)=\frac d h,\qquad v(S)=-\frac{d+2}{h}.
\]

The nonvanishing derivative of (15) gives

\[
 v(D)=-\frac{d+4}{h},\qquad
 v(q)=v(PD)=-\frac4h<0.                                 \tag{16}
\]

The polynomial source coordinate `q` therefore has a pole.  Finally,

\[
 2+1+h e_0=N                                            \tag{17}
\]

is the full local degree sum over `P=0`.

### Exhaustion and intrinsic ordering

If `P` and the reduced discriminant equation are both nonzero, every inverse
root is simple and the reconstruction formulas

\[
 t=D^{-1},\qquad x=S/D,\qquad q=PD
\]

together with the polynomial formulas for `y` and `z` are regular.  Thus
every divisorial boundary image is contained in `Z_Delta union Z_0`, where
`Z_0=V(P)`.  Equations (9) and (17) exhaust the generic degree over both
divisors.  Purity of the complement in the affine finite normalization, as
in the canonical boundary-exhaustion theorem, excludes a residual
higher-codimension component.

Hence (12) and the discriminant row are the complete canonical boundary
list.  The two target vertices are intrinsically ordered.  The vertex
`Z_Delta` receives exactly one boundary prime with `(e,f)=(2,1)`.  The
profile over `Z_0` is:

\[
 \begin{cases}
  \text{one }(1,1),&N=4,\\
  \text{two }(1,1),&N=5,\\
  h\text{ copies of }(d/h,1),&N\ge6.
 \end{cases}                                            \tag{18}
\]

This never equals the discriminant profile.  In particular, when `d=4`
the second vertex receives two index-two primes, not one.  Ramification and
residue degrees and the number of primes are part of the map-intrinsic
normalization package, so the ordering does not depend on the displayed
coordinate `P`.

For comparison, the weighted clean-locus theorem gives the same ordered
pair notation: its `Z_Delta` receives one `(2,1)` boundary prime and its
second vertex receives `N-3` unramified `(1,1)` boundary primes.  Already
the second-vertex ledgers separate the two families for every `N>=6`.

## 2. Normalization of the intrinsic ramified stratum

Delete the intrinsically selected second vertex.  For the quadratic family
this means `P!=0`.  Let

\[
 R_\Delta=
 k[P^{\pm1},B(r),C(r)]
 \subset k[P^{\pm1},r^{\pm1}].                          \tag{19}
\]

The equation

\[
 G_P'(r)-g_1B r=0                                      \tag{20}
\]

has leading coefficient `Ng_NP^N`, a unit in `k[P^{\pm1}]`.  After division
by that unit it is monic in `r`, so `r` is integral over `R_Delta`.
Moreover,

\[
 r^{-1}
 =B-\frac{G_P'(r)-g_1}{g_1r},                           \tag{21}
\]

and the quotient on the right is polynomial in `r`.  Thus `r^{-1}` belongs
to `R_Delta[r]`.  It follows that

\[
 k[P^{\pm1},r^{\pm1}]
\]

is finite over `R_Delta`.  It is normal, and (7) says that the two rings have
the same fraction field.  Therefore

\[
 \boxed{
 \operatorname{Norm}(Z_\Delta\setminus Z_0)
 \simeq\operatorname{Spec}k[P^{\pm1},r^{\pm1}]
 =\mathbb G_m^2.
 }                                                       \tag{22}
\]

This calculation also proves that both missing points `r=0,infinity` remain
missing in the canonical affine stratum: neither is filled by passing from
the fixed-`P` pencil to the full three-dimensional normalization.

For the weighted clean locus, the intrinsic-boundary theorem gives

\[
 \boxed{
 \operatorname{Norm}(Z_\Delta\setminus Z_0)
 \simeq
 \operatorname{Spec}k[W,\xi^{\pm1}]
 =\mathbb A^1\times\mathbb G_m.
 }                                                       \tag{23}
\]

Both (22) and (23) were obtained only after selecting and deleting canonical
vertices of the Zariski--Main package.  They are therefore intrinsic finite
strata, not presentation-dependent slices.

## 3. Stable unit-rank obstruction

The construction-independent
[stable normalization functoriality theorem](STABLE_NORMALIZATION_FUNCTORIALITY.md)
transports the ordered boundary package under polynomial left--right
equivalence.  Under stabilization it replaces each stratum by its product
with `A^s`.

For every `s>=0`,

\[
 \begin{aligned}
 \mathcal O(\mathbb A^1\times\mathbb G_m\times\mathbb A^s)^\times
   &=k^\times\xi^{\mathbb Z},\\
 \mathcal O(\mathbb G_m^2\times\mathbb A^s)^\times
   &=k^\times P^{\mathbb Z}r^{\mathbb Z}.               \tag{24}
 \end{aligned}
\]

Their unit groups modulo `k^\times` have ranks one and two.  An isomorphism of
stabilized intrinsic strata would preserve the unit group, which is
impossible.  This proves the stable separation theorem.

## 4. Scope of the word “components”

The theorem proves disjointness of the two loci in the quotient by stable
polynomial left--right equivalence: no stable orbit meets both.  Calling
them distinct irreducible or connected components of a moduli space requires
a separately constructed moduli space and a proof that the unit-rank strata
are open-and-closed there.  No such additional moduli-topology claim is
needed for, or made by, the stable inequivalence theorem.

## Exact regression

Run

```bash
.venv/bin/python scripts/verify_quadratic_weighted_stable_separation.py
```

The checker verifies the inverse derivative, the two affine `P=0` branches,
the Newton polygons and degree sums in a broad degree range, the high-branch
source pole, and the Laurent recovery identities used in (20)--(22).
