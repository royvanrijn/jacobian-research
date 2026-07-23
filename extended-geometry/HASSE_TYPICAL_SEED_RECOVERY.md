# Minimal p-typical Hasse seed recovery

The ordinary Hessian portrait is not coefficient-faithful in bad
characteristic.  This is not a boundary pathology: Frobenius directions can
move inside the Hessian-clean locus while leaving the projective ordinary
Hessian unchanged.  The degree-eight collision over `F_5` in
[`D1_F2_COEFFICIENT_SPACE_PROOF.md`](D1_F2_COEFFICIENT_SPACE_PROOF.md)
is an explicit example.

There is nevertheless a small and sharp coefficient-level repair.  The
ordinary second derivative is replaced by the Hasse derivatives in the
`p`-power orders

```text
1, p, p^2, ..., p^e,       e=floor(log_p(N)).
```

Only `e+1` polynomial channels are needed.  This repairs coefficient
recovery, but it does not by itself construct a map-intrinsic
positive-characteristic stable invariant.  Section 6 records that boundary
precisely.

## 1. Hasse derivatives

Let `k` be a field of characteristic `p>0`.  The `m`-th Hasse derivative is
the `k`-linear operator determined by

```text
D^[m](W^j) = binom(j,m) W^(j-m),                       (1.1)
```

with the convention that the right side is zero for `j<m`.  Equivalently,

```text
H(W+T) = sum_{m>=0} D^[m]H(W) T^m.                    (1.2)
```

In particular, `D^[1]` is the ordinary first derivative and

```text
H'' = 2 D^[2]H.                                       (1.3)
```

Formula (1.3) also explains why an ordinary Hessian is especially lossy in
small characteristic: it is only one Hasse channel, and in characteristic
two it vanishes identically.

More precisely,

```text
ker(d^2/dW^2) = k[W^p] direct_sum W k[W^p].           (1.4)
```

Indeed, the coefficient multiplying `W^(j-2)` is `j(j-1)`, which vanishes
exactly when `j` is congruent to zero or one modulo `p`.  Thus the ordinary
Hessian forgets two entire Frobenius residue classes before any projective
quotient is taken.

## 2. The reconstruction theorem

Fix `N>=1`, put

```text
e = floor(log_p(N)),
T_(p,N)(H) = (D^[1]H, D^[p]H, ..., D^[p^e]H).         (2.1)
```

### Theorem

The linear map

```text
T_(p,N): k[W]_(<=N) -> direct_sum_(i=0)^e k[W]_(<=N-p^i)
```

has kernel exactly `k`.  Consequently it is injective modulo constants.  On
any normalized seed space on which `H(0)=0`, the tuple (2.1) recovers `H`
exactly.

More explicitly, write

```text
H(W)=sum_(j=0)^N a_j W^j,
j=sum_(r>=0) j_r p^r,       0<=j_r<p.                 (2.2)
```

For every `j>0`, choose any index `i(j)` with `j_(i(j))!=0`.  Then

```text
a_j =
  j_(i(j))^(-1)
  [W^(j-p^i(j))] D^[p^i(j)]H.                         (2.3)
```

Thus (2.3), together with a supplied value of `a_0`, is a linear left
inverse to (2.1).  Choosing the least nonzero base-`p` digit gives a
canonical reconstruction rule.

### Proof

The characteristic-`p` identity

```text
(1+T)^j
  = product_(r>=0) (1+T^(p^r))^(j_r)                  (2.4)
```

is the special case of Lucas' theorem needed here.  Taking the coefficient
of `T^(p^i)` gives

```text
binom(j,p^i) = j_i in k.                              (2.5)
```

In `D^[p^i]H`, the coefficient of `W^(j-p^i)` can only come from the
monomial `a_jW^j`; no other input exponent has the same shifted exponent.
If `j_i!=0`, equations (1.1) and (2.5) therefore give

```text
[W^(j-p^i)]D^[p^i]H = j_i a_j.                        (2.6)
```

The element `j_i` lies in `{1,...,p-1}`, so it is invertible in `k`.
Every positive `j` has a nonzero base-`p` digit, and (2.3) follows.
All positive-degree coefficients are consequently recovered.  Constants
are killed by every positive Hasse derivative, proving that the kernel is
exactly `k`.

The same proof works over every commutative `F_p`-algebra: each nonzero
base-`p` digit is still a unit.  In that form, (2.1) is a split injection of
finite free modules after quotienting the source by constants, so the result
is scheme-theoretic and stable under base change.

## 3. Sharp minimality among Hasse orders

The channel list in (2.1) is not merely sufficient.

### Minimality theorem

Let `S` be any set of positive Hasse-derivative orders at most `N`.  If

```text
H |-> (D^[m]H)_(m in S)                               (3.1)
```

is injective on `k[W]_(<=N)` modulo constants, then

```text
{1,p,p^2,...,p^e} is a subset of S.                   (3.2)
```

Hence every such family has at least `e+1` channels, and the p-typical
family is the unique minimum-cardinality family.

### Proof

Fix `p^i<=N`.  Lucas' theorem applied to the base-`p` expansion of `p^i`
gives

```text
binom(p^i,m) = 0 in k       for 0<m<=N, m!=p^i.       (3.3)
```

For `m>p^i`, the derivative vanishes by degree.  Therefore

```text
D^[m](W^(p^i))=0            for every positive m!=p^i,
D^[p^i](W^(p^i))=1.                                  (3.4)
```

If `p^i` is omitted from `S`, the nonconstant polynomial `W^(p^i)` lies in
the kernel of (3.1), contradicting injectivity modulo constants.

This minimality statement is deliberately about raw Hasse-derivative orders
on the full space `k[W]_(<=N)/k`.  It does not assert that the same number is
always minimal after restricting to a smaller normalized-seed affine slice:
relations such as `H(1)=0` can make a coefficient dependent on the others.
The p-typical family remains a uniform sufficient family on every such
slice.

## 4. Degree eight over F_5

For `p=5` and `N=8`, the theorem needs only

```text
D^[1]H and D^[5]H.                                   (4.1)
```

Recall the two normalized, exact-degree, primitive-squarefree, and
Hessian-squarefree seeds

```text
H_1 = W^2+4W^3+W^6+3W^7+W^8,
H_2 = 2W^2+3W^3+4W^5+3W^6+W^7+2W^8.                 (4.2)
```

Their ordinary Hessians are different scalar representatives of the same
projective polynomial (the displayed tuples list coefficients from low to
high):

```text
H_1'' = (2,4,0,0,0,1,1),
H_2'' = (4,3,0,0,0,2,2) = 2H_1'',                   (4.3)
[H_1'']=[H_2'']=[1,2,0,0,0,3,3].
```

In fact, with `L_5=W^5-W^6`,

```text
H_2 = 2H_1-L_5,       L_5''=0.                        (4.4)
```

This identifies the collision direction exactly.  The ordinary marked
Hessian portrait collides, but the missing Frobenius coefficient is visible
immediately in the fifth Hasse channel:

```text
D^[5]H_1 = W+3W^2+W^3,
D^[5]H_2 = 4+3W+W^2+2W^3.                            (4.5)
```

Together with `D^[1]`, formula (2.3) recovers both seeds exactly.  At the
level of individual monomials, `W^5` is invisible to every positive Hasse
order other than five, which is the local reason that the fifth channel
cannot be replaced.

The collision (4.2) involves projective ordinary Hessians, whereas the
reconstruction theorem uses actual polynomial channels.  This distinction
is essential for the geometric application.

The same mechanism repairs the full positive-dimensional family from the
bad-characteristic audit.  For

```text
L_p=W^p-W^(p+1),       H_c=cG+(1-c)L_p,
```

one has `H_c''=cG''`, so `[H_c'']` is constant for `c!=0`, while

```text
D^[p]L_p=1-W.                                         (4.6)
```

If two members of this family have the same complete p-typical tuple, the
reconstruction theorem makes their difference constant; their common
normalization at zero then makes them equal.  Hence this
positive-dimensional marked Hessian fiber disappears at the coefficient
level.

## 5. Covariance and the natural weighted package

Hasse derivatives have clean affine covariance.  For `a!=0`,

```text
D^[m](H(aW+b)) = a^m (D^[m]H)(aW+b).                 (5.1)
```

The p-typical tuple therefore transforms with channel weights

```text
1,p,p^2,...,p^e.                                      (5.2)
```

Multiplying `H` by one scalar multiplies every channel by that same scalar.
Consequently, retaining the tuple in one direct sum preserves meaningful
relative channel scales.  A common projectivization of that direct sum loses
only the common scalar; on a seed normalized by `H'(1)=-1`, that scalar can
be restored from the first channel.

By contrast, replacing every nonzero channel independently by its divisor,
or independently projectivizing the channels, forgets their relative
scalars.  Some channels may also vanish on special seed strata.  Either
operation requires additional analysis before coefficient recovery descends
through a quotient.

## 6. A map-intrinsic no-go theorem for the existing construction

The remaining gap cannot be closed merely by applying a more elaborate
functor to the existing weighted map.  The map itself forgets some of the
coefficients recovered by the p-typical tuple.

The weighted-map formula is built from the ordinary seed

```text
P(W)=H'(W)
```

and from expressions obtained from `P` and `P'`.  Consequently it factors
through the ordinary derivative `H |-> H'` in every characteristic where
the reduced formula is defined.

### No-go theorem

Let `k` have characteristic `p>0`, and suppose `H` is normalized by

```text
H(0)=H'(0)=H(1)=0,       H'(1)=-1.                    (6.1)
```

If `K in k[W^p]` satisfies `K(0)=K(1)=0`, then every member

```text
H_c=H+cK.                                             (6.2)
```

has the same ordinary derivative and the same endpoint normalization.
Whenever the reduced weighted formula for their common derivative is
polynomial, the existing construction assigns them the **same polynomial
map**.  Therefore
no invariant reconstructed functorially from that map—not even the entire
Zariski--Main normalization package with all of its finite strata—can
recover `H_c`.

On the other hand, the p-typical tuple separates distinct members modulo
constants.  Thus the coefficient repair in Section 2 cannot be promoted to
a map-intrinsic invariant of the existing construction.  The construction
must first be enriched so that its polynomial coordinates depend on the
missing Hasse channels.

### Proof

Every exponent occurring in `K` is divisible by `p`, so `K'=0`.  Hence

```text
H_c'=H'.                                               (6.3)
```

The endpoint conditions on `K` preserve (6.1).  The formulas defining the
weighted map use `H'`, its ordinary derivatives, and the fixed construction
constants; all are independent of `c`.  The resulting polynomial maps are
therefore literally equal, not merely stably left--right equivalent.

If two p-typical tuples were equal, Section 2 would make
`H_c-H_d=(c-d)K` constant.  Its value at zero is zero, so it would vanish.
For nonzero `K`, this forces `c=d`.

There is also a no-go statement already at the discriminant-normalization
level.  Put `f=H'` and

```text
R(W)=(f(W)(f(W)+1))^p.                                (6.4)
```

Then `R(0)=R(1)=0`, `R'=0`, and the repeated-root parameterizations satisfy

```text
nu_(H+cR)(W)
 = (f(W), Wf(W)-H(W)-cR(W))
 = tau_c(nu_H(W)),                                   (6.5)

tau_c(s,t)=(s,t-c(s(s+1))^p).
```

The map `tau_c` is a triangular target automorphism and the normalization
coordinate, including the marks `0`, `1`, and infinity, is unchanged.
Thus the complete marked finite normalization morphism alone also cannot
supply the missing channel.

### A clean collision in degree twelve over F_5

The obstruction occurs inside the same clean conditions used by the
ordinary Hessian audit.  Over `F_5`, put

```text
H_0 =
  3W^2+2W^3+4W^4+2W^5+2W^6+W^7
  +2W^8+W^9+3W^10+3W^11+2W^12,
K = W^10-W^5.                                        (6.6)
```

For

```text
c in F_5,
```

the five seeds `H_c=H_0+cK` all have:

* exact degree twelve;
* the normalization (6.1);
* squarefree `H_c/W^2` of degree ten;
* the same squarefree ordinary Hessian of degree ten;
* the same ordinary derivative, hence the same weighted polynomial map.

Their common Hessian coefficient vector, from low to high, is

```text
(1,2,3,0,0,2,2,2,0,0,4).                             (6.7)
```

But

```text
D^[5]K = 4+2W^5 != 0,                                (6.8)
```

so the minimal p-typical tuple distinguishes all five seeds.  This is a
clean, same-degree obstruction to map-intrinsic seed recovery, stronger than
the degree-eight projective-Hessian collision.

## 7. A dimension-preserving Frobenius-enriched Keller lift

Although the existing map cannot remember `K`, there is a direct polynomial
enrichment that makes it do so.  This closes the first item in the
replacement program without adding variables.

Write the weighted source expressions as

```text
u=1+xy,
gamma=1+a xy+b x^2z,
W=u gamma,
C=x gamma.                                            (7.1)
```

Suppose the ordinary-derivative construction gives a polynomial Keller map

```text
G_0=(A_0,B,C)
```

whose plane-core inverse polynomial is `H_0`, and let
`K in k[W^p]` satisfy `K(0)=K(1)=0`.

### Enrichment theorem

The quotient

```text
Delta_K(x,y,z)=K(W)/C^2                               (7.2)
```

is a polynomial.  Moreover

```text
G_K=(A_0-Delta_K/c,B,C)                               (7.3)
```

has the same constant Jacobian determinant as `G_0`, and its inverse pencil
on `C!=0` is

```text
H_0(W)+K(W)-BCW+cAC^2=0.                              (7.4)
```

Thus the correction retains the ambient dimension, preserves the Keller
property, and makes the map depend on the full Frobenius primitive.

### Proof

Write `K(W)=bar K(W^p)`.  The endpoint conditions imply

```text
bar K(T)=T(T-1)Q(T).
```

Since `W-1=u gamma-1` is divisible by `x`,

```text
K(W)
 = W^p(W-1)^p Q(W^p)
 = u^p gamma^p x^p ((W-1)/x)^p Q(W^p).
```

After division by `C^2=x^2gamma^2`, the remaining factors
`x^(p-2)` and `gamma^(p-2)` have nonnegative exponents.  This proves
polynomiality for every `p>=2`.

The target expression `cAC^2` changes by exactly `-K(W)`.  Hence the
weighted suspension square now has plane core

```text
(W,gamma) |->
  (H_0'(W)+c gamma,
   W(H_0'(W)+c gamma)-H_0(W)-K(W)).
```

Because `K'=0`, its Jacobian is the same as the original plane core.
The determinant cancellation in the weighted suspension is unchanged,
proving both the Keller assertion and (7.4).

### The enriched degree-twelve family

For (6.6), the ordinary-derivative map corresponds to the representative
whose `W^10-W^5` coefficient is zero.  The five corrected first coordinates
are

```text
A_c=A_ordinary-(3+c)(W^10-W^5)/C^2,       c in F_5.   (7.5)
```

The correction is a degree-42 polynomial with 60 monomials.  Exact expansion
gives five distinct polynomial maps with coordinate degrees `(52,51,4)`,
Jacobian determinant one, and inverse pencils belonging to the five
polynomials `H_c`.

The existing discriminant decorations see the full enrichment in this
explicit family.

First, the family-specific intrinsic-boundary argument remains valid in
characteristic five.  The discriminant branch is generically a tame double
root; the zero cluster is also tame of degree two; `H_c/W^2` and `H_c''`
are squarefree; and the normalization coordinates have coprime pole orders
`11` and `12`.  The enriched correction preserves the plane-core
reconstruction identities.  The degree sums over the discriminant and zero
boundaries are therefore still `2+10=12` and `2+1+9=12`; as in the
characteristic-zero proof, they exclude hidden boundary primes.
Consequently the boundary pair, normalization line, and root-one affine
sheet are intrinsic for these five maps.  The three marks `0`, infinity, and
`1` kill every automorphism of the normalization line.

For a quick visible distinction, divide the two root-one equal-image
equations by `u-1` and take their gcd.  Its degree for
`c=0,1,2,3,4` is respectively

```text
0,0,1,0,0.                                            (7.6)
```

At `c=2` the gcd is `u+2`, so the marked point `1` shares its discriminant
image with `u=3`.  The Hessian values there are `1` and `3`, and the tangent
directions are `(1,1)` and `(1,3)`, proving that this is a transverse node.
The other four members have no geometric root-one self-pair.

The complete equal-image schemes separate the other four maps as well.
After removing the common diagonal factor, compute the reduced lexicographic
Groebner basis in normalization coordinates `(r,u)`.  Its final monic
eliminant has degree `110`.  Listing its coefficients of
`1,u,u^2,u^3` gives

| `c` | low-coefficient fingerprint |
|---:|---:|
| `0` | `(4,1,4,0)` |
| `1` | `(3,0,2,2)` |
| `2` | `(3,3,4,4)` |
| `3` | `(2,1,2,3)` |
| `4` | `(2,3,4,2)` |

These are five distinct presentations of the intrinsic normalization
self-fiber product.  Since every isomorphism of the three-marked
normalization line is the identity, no two schemes can be identified.
Stable normalization functoriality preserves the finite stratum, its fiber
product, and the intrinsic marks in arbitrary characteristic.  Therefore:

> **Explicit stable-separation theorem.**  The five Frobenius-enriched
> degree-twelve Keller maps (7.5) over `F_5` are pairwise stably polynomial
> left--right inequivalent.

### Tame all-degree intrinsic faithfulness

The same mechanism is not restricted to the example.  Let `p` be odd and
consider the Frobenius-enriched maps (7.3) on the open where:

1. `H` has exact degree `N`, an exact double root at zero, and all other
   primitive roots are simple;
2. `H''` is squarefree of exact degree `N-2` with at least two support
   points;
3. the weighted reconstruction is boundary-clean.

These hypotheses exclude `p|N(N-1)` and every wild degree-two branch.

> **Tame intrinsic-faithfulness theorem.**  On this open, the enriched
> polynomial map determines the normalized seed `H` exactly under stable
> polynomial left--right equivalence.

Indeed, the intrinsic boundary profile retains every extra-root prime over
the second target boundary and its finite residue factor.  In the
normalization coordinate, the product of those factors is

```text
H(W)/(W^2(W-1)),                                      (7.7)
```

up to one scalar.  This collection is unordered, as it must be, but its
product is canonical.  The full Fitting divisor has at least two support
points, so the Laurent-unit argument forces every stabilized
normalization-line change to be affine.  The intrinsic zero and infinity
marks reduce it to a scaling, and the intrinsic affine root-one sheet forces
that scaling to be one.  Finally `H'(1)=-1` fixes the scalar in (7.7).

All ingredients are functorial parts of the Zariski--Main boundary package:
the ordered target boundaries, finite edge maps and residue factors, full
Fitting divisor, and affine sheet.  Tameness makes the
characteristic-zero boundary proof valid verbatim.  Hence a stable
left--right equivalence of two enriched maps forces equality of their
normalized seeds.

This proof deliberately uses the stronger full boundary invariant rather
than trying to turn each Hasse channel separately into a divisor.  The
p-typical theorem reconstructs and parametrizes the input coefficients; the
enriched map and its complete boundary recover those coefficients stably.

### Characteristic two: exact obstruction and residual data

Characteristic two is not just a wild version of the preceding proof.  If

```text
H(W)=sum_m a_(2m)W^(2m)+sum_m a_(2m+1)W^(2m+1),
```

put

```text
A(T)=sum_m a_(2m+1)T^m,
B(T)=sum_m a_(2m)T^m.                                (7.8)
```

Then

```text
H'(W)=A(W^2),
WH'(W)-H(W)=B(W^2),
H(W)=W A(W^2)+B(W^2).                                (7.9)
```

Thus the critical parameterization factors through the Frobenius coordinate
`T=W^2`.  It is not the normalization map from the `W`-line; the reduced
discriminant normalization must first be rebuilt on the `T`-line.
Nevertheless the two parity polynomials `(A,B)` still recover `H` exactly,
so the coefficient data themselves have not disappeared.

The obstruction occurs earlier in the weighted suspension.  Its defining
parameter is

```text
a=-(1+kappa)/(2+kappa),       kappa=H''(1)/c.          (7.10)
```

In characteristic two, `H''=0` for every polynomial.  Hence
`kappa=0` and the denominator in (7.10) is zero identically.  The existing
weighted construction has no characteristic-two reduction to enrich.
A new suspension based on the Frobenius-compressed `(A(T),B(T))` core is
required.

This failure persists under the natural scalar reparameterizations of the
same suspension.  Allow

```text
W=(1+dxy)gamma,
s=H'(W)+cR(gamma),
C=xgamma,                                             (7.11)
```

with `R(0)=0`, `R(1)=1`.  Divisibility of the second target coordinate by
`C` requires the first condition, while the root-one normalization requires
the second.  The plane-core determinant is

```text
-c^2 R(gamma)R'(gamma).
```

After the two vertical Jacobians are included, constancy of the threefold
Jacobian requires

```text
R(gamma)R'(gamma)/gamma = nonzero constant.           (7.12)
```

Writing `R=gamma S`, equation (7.12) says that the product
`S(S+gamma S')` is a nonzero constant.  Both factors must therefore be
constant.  The endpoint normalization then forces `R=gamma`.

For `R=gamma`, the coefficient of `xy` in the numerator that must be
divisible by `x^2` is

```text
c((2+kappa)a+(1+kappa)d)=cd                          (7.13)
```

in characteristic two.  Polynomiality forces `d=0`.  But the source
vertical Jacobian in (7.11) is proportional to `d`, so `d=0` makes the
suspension degenerate.  Thus:

> **Characteristic-two ansatz no-go theorem.**  No separable scalar
> reparameterization (7.11) simultaneously gives polynomial weighted
> coordinates and a nonzero constant Jacobian.  A characteristic-two
> continuation must change the plane core or the suspension architecture.

### Characteristic two: a weight-redistributed suspension

There is a smaller change that escapes the no-go theorem: redistribute the
two source weights rather than reparameterizing the plane scalar.  Assume
`char(k)=2` and

```text
H(0)=H'(0)=H(1)=0,       H'(1)=1.
```

Put

```text
u=1+x^2y,        gamma=1+xz,
W=u gamma,       C=x gamma,
S=H'(W)+gamma,
T=W S-H(W).                                         (7.14)
```

> **Characteristic-two suspension theorem.**  Both quotients in
>
> ```text
> G_H=(T/C^2,S/C,C)                                  (7.15)
> ```
>
> are polynomials, and `det J(G_H)=1`.  On `C!=0` its inverse pencil is
>
> ```text
> H(W)-BCW+AC^2=0,       (A,B,C)=G_H.                 (7.16)
> ```
>
> Thus normalized characteristic-two seeds admit a dimension-preserving
> polynomial Keller realization.  Unlike the ordinary-derivative
> construction, (7.15) depends on the complete seed `H`.

To prove polynomiality, use the parity split (7.8).  Since
`H'=A(W^2)`,

```text
S=A(u^2 gamma^2)+gamma.
```

The condition `H'(0)=0` makes `S` divisible by `gamma`, while its value at
`x=0` is `H'(1)+1=0`; hence `S/(x gamma)` is polynomial.  Also

```text
T=u gamma^2-B(u^2 gamma^2).                           (7.17)
```

Every term is divisible by `gamma^2`, because `B(0)=0`.  After that factor
is removed, the result is a polynomial in `u` and `gamma^2`.  Its value at
`(u,gamma)=(1,1)` is zero: (7.9) and the endpoint conditions give
`B(1)=A(1)=1`.  It therefore lies in

```text
(u-1,gamma^2-1)=(x^2y,x^2z^2),
```

which proves divisibility by `x^2 gamma^2=C^2`.

For the determinant, the three factors in the suspension square are

```text
det d(W,gamma,C)/d(x,y,z) = x^3 gamma^2,
det d(T,S)/d(W,gamma)     = gamma,
det d(T/C^2,S/C,C)/d(T,S,C) = C^(-3).                (7.18)
```

Their product is one because `C=x gamma`.  Equation (7.16) is the defining
identity for `T`.  Notice why this construction is outside (7.11):
`u-1` begins in order `x^2`, so the polynomiality obstruction vanishes,
whereas `gamma-1` begins in order `x`; their vertical Jacobian still has
the required order `x^3`.

The weight redistribution is actually characteristic-free on exactly the
singular locus of the old parameter.  For a normalized seed in any
characteristic, put `kappa=H''(1)`.  The derivative with respect to `gamma`
of `T/gamma^2` at `(u,gamma)=(1,1)` is

```text
2+kappa.                                             (7.18a)
```

Hence (7.14)--(7.18) remain valid whenever `2+kappa=0`: the numerator lies
in `(u-1,(gamma-1)^2)` and is divisible by `x^2`.  When
`2+kappa!=0`, the original chart uses

```text
u=1+xy,
gamma=1-((1+kappa)/(2+kappa))xy+x^2z.
```

Thus the two source charts cover every normalized seed for which the plane
core is defined.  The integral quartic

```text
H=2W^2-3W^3+W^4
```

has `H''(1)=-2`; the checker verifies the redistributed determinant-one map
for its reductions in characteristics three, five, and seven.

If `H` has exact degree `N`, then (7.15) has geometric degree `N`.  Indeed,
for a target point `(a,b,c)` with `c!=0`, every root of

```text
H(W)-bcW+ac^2=0                                      (7.19)
```

for which `gamma=bc-H'(W)` is nonzero gives exactly one source point:

```text
x=c/gamma,   u=W/gamma,
y=(u-1)/x^2, z=(gamma-1)/x.                          (7.20)
```

Conversely these formulas recover every source point.  The derivative of
(7.19) is `H'(W)-bc=-gamma`, so the generic degree-`N` fiber is separable
despite the Frobenius factorization of the critical-value curve.

For example, the normalized cubic `H=W^2+W^3` gives the explicit
characteristic-two Keller map

```text
(x,y,z) |->
(
  y+x^2y^2,
  z+x^3y^2+x^4y^2z,
  x+x^2z
),                                                    (7.21)
```

with Jacobian one and geometric degree three.

### Characteristic two: wild-clean intrinsic faithfulness

The reduced discriminant does not retain the `W`-parameter as its
normalization coordinate, but the full boundary package retains the missing
purely inseparable edge.  For a new variable `R`, specialize the inverse
pencil to the compressed discriminant parameter:

```text
s=A(T),       t=B(T).
```

Then

```text
H(R)-A(T)R+B(T)
 = R(A(R^2)-A(T))+B(R^2)-B(T),                       (7.22)
```

which is divisible by `R^2-T`.  Generically this is the complete
double-root cluster.  Thus the boundary edge over the reduced discriminant
normalization is the finite radicial morphism

```text
A^1_W -> A^1_T,       T=W^2,                         (7.23)
```

of residue degree two.  The other `N-2` roots are regular source sheets.

At the second target boundary `C=0`, write `W=CR` and
`H(W)=h_2W^2+O(W^3)`.  The zero cluster has special equation

```text
h_2R^2-BR+A=0.                                       (7.24)
```

Its derivative is `B`, so it is generically separable of residue degree
two.  This intrinsically orders the two target boundaries: one carries a
radicial degree-two edge and the other a separable degree-two affine
cluster.  The root `W=1` gives the unique additional residue-degree-one
affine sheet.  Indeed, the reconstruction formulas (7.20) show that a
simple root `rho` can remain affine only if first
`H'(rho)=rho` and then `H'(rho)=1`, hence `rho=1`.

The birationality of the compressed map is automatic on the primitive
squarefree locus.  Suppose instead that
`k(A(T),B(T))` is a proper subfield of `k(T)`.  The polynomial form of
Luroth's theorem gives a polynomial `Q` of degree `d>1` and polynomials
`A_0,B_0` such that

```text
A=A_0 o Q,       B=B_0 o Q.                           (7.24a)
```

At `q_0=Q(0)`, both `A_0(q_0)` and `B_0(q_0)` vanish.  Every nonzero root
`tau` of `Q(T)-q_0`, followed by a square root `rho^2=tau`, then satisfies

```text
H(rho)=H'(rho)=0.
```

If no such nonzero `tau` exists, then
`Q(T)-q_0=cT^d`; but every polynomial in `Q` has zero linear `T`
coefficient, contradicting the nonzero `T` coefficient `h_2` of `B(T)`
coming from the exact double zero.  Thus a proper compressed subfield forces
a repeated nonzero primitive root.

Call a normalized seed **wild-clean** if:

1. it has exact degree `N>=4`, an exact double root at zero, and all its
   other roots are simple;
2. regular reconstruction has no further codimension-one failure.

> **Characteristic-two wild-clean theorem.**  On this locus, the full
> intrinsic boundary package of (7.15) determines the normalized seed `H`
> exactly under stable polynomial left--right equivalence.

Over the discriminant boundary, the radicial cluster contributes degree two
and the regular sheets contribute `N-2`.  Over `C=0`, (7.24), the root-one
sheet, and the extra-root boundary primes contribute `2+1+(N-3)=N`.
These degree sums, followed by the normal Hartogs argument, exclude hidden
height-one and higher-codimension boundary components exactly as in the
tame proof.

The radicial edge retains the affine `W`-line and its zero-cluster and
root-one marks.  After stabilization these are the two principal marked
ideals `(W)` and `(W-1)`.  If an automorphism carries them to the
corresponding two ideals, the polynomial-ring unit theorem gives

```text
W_new=aW,       W_new-1=b(W-1)
```

for constants `a,b`; comparison forces `a=b=1`.  Thus the stabilized
boundary package fixes `W` itself without any Hessian-support hypothesis.
The finite residue factors on the extra-root primes have
product

```text
H(W)/(W^2(W-1))                                      (7.25)
```

up to one scalar.  The fixed coordinate and `H'(1)=1` recover that scalar
and hence `H`.  This proves the theorem.  The argument uses the radicial
edge itself in place of the identically zero ordinary-Hessian Fitting
divisor.

The wild-clean locus is nonempty in every degree `N>=5`.  Over `F_2`, use

```text
H_N=W^2+W^N                         for odd N,
H_N=W^2+W^3+W^4+W^N                for even N.        (7.26)
```

For odd `N`, the compressed pair is
`(A(T),B(T))=(T^((N-1)/2),T)`; for even `N`, it is
`(T,T+T^2+T^(N/2))`.  Hence the compressed map is birational.  Moreover
`H_N/W^2` is squarefree: its derivative is `W^(N-3)` in odd degree and
`1` in even degree.  The checker verifies (7.22), the exact degrees, and
these squarefreeness claims through degree sixteen; the formulas prove them
uniformly.

Degree three has no residual moduli: the normalization conditions force

```text
H=W^2+W^3.
```

Thus (7.21) also gives stable seed-faithfulness in degree three
tautologically.  Over a characteristic-two field the exact-degree-four
normalized seeds are

```text
H_lambda=(1+lambda)W^2+W^3+lambda W^4,
lambda!=0,1.                                         (7.27)
```

For `lambda!=0,1`, the primitive quotient is squarefree, the compressed
pair has first coordinate `A(T)=T`, and the only extra-root residue factor
is

```text
lambda W+(1+lambda).
```

Consequently this entire family is wild-clean and (7.25) recovers
`lambda`.  Characteristic-two stable seed-faithfulness is therefore
complete on the clean locus in every degree `N>=3`.

### Unified full-edge clean theorem

The preceding marked-edge argument also removes the Hessian degree and
support restrictions from the odd-characteristic theorem.  Let `k` have
arbitrary positive characteristic, let `H` have an exact double zero and
only simple nonzero roots, and use the original or redistributed source
chart according as `2+H''(1)` is nonzero or zero.  Assume that the critical
boundary parameterization is generically birational onto its reduced image
and that regular reconstruction is boundary-clean.

> **Full-edge clean theorem.**  The complete intrinsic boundary package
> determines `H` exactly under stable polynomial left--right equivalence.
> No exact-degree or two-support hypothesis on `H''` is required.

In odd characteristic the critical boundary edge itself is the `W`-line;
in characteristic two it is the radicial edge (7.23).  In both cases its
zero-cluster and root-one intersections give the two intrinsic principal
ideals `(W)` and `(W-1)`, which fix `W` after arbitrary polynomial
stabilization.  The second-boundary residue factors again multiply to
(7.25), and endpoint normalization fixes their common scalar.  The same
degree sums prove completeness of the two boundary targets.  The
Hessian/Fitting divisor is therefore a useful smaller invariant on the tame
open, but it is not needed once the full finite edge data are retained.

The birationality hypothesis is automatic in odd characteristic.  Put

```text
u=H'(W),       v=WH'(W)-H(W).
```

In the one-variable function field,

```text
dv=W du.                                             (7.28)
```

If `H''!=0`, the ratio `dv/du` therefore places `W` in `k(u,v)`, proving
`k(u,v)=k(W)`.  For an exact double zero
`H=h_2W^2+O(W^3)`, odd characteristic gives
`H''(0)=2h_2!=0`.  In characteristic two, compressed birationality on the
clean locus was proved above.  Consequently the full-edge clean theorem
applies to every normalized boundary-clean seed with an exact double zero
and simple nonzero roots, in every positive characteristic.

### Collision extension and the smaller exact quotient

The simplicity hypothesis can also be removed.  Let `rho!=0,1` be a root
of multiplicity `m`, and write

```text
H(rho+delta)=h_m delta^m+O(delta^(m+1)),       h_m!=0.
```

At the generic point of the second target boundary `C=0`, the inverse
pencil has local equation

```text
h_m delta^m-B rho C-B C delta+A C^2
  +O(delta^(m+1))=0.                                (7.29)
```

Its Weierstrass polynomial is Eisenstein at `C`: the constant coefficient
has valuation one and the leading coefficient is a unit.  The normalized
boundary therefore has one prime of ramification index `m` centered at
`W=rho`.  This remains true when `p|m`; the edge is then wild or radicial,
but its valuation and residue center are unchanged.  Over a nonclosed base
field, the same statement retains the irreducible residue factor and its
multiplicity.

The zero root behaves differently and remains intrinsic.  Since its
multiplicity is exactly two, `W=CR` gives the finite affine cluster (7.24).
The root `1` remains the unique other affine sheet.  All remaining roots,
with their full multiplicities, are precisely the normalized boundary
primes over `C=0`.  Their degree sum is

```text
2+1+sum_(rho!=0,1) mult_rho(H)=N.                    (7.30)
```

Thus no hidden boundary prime is introduced at a collision.
Although a repeated extra root makes this target boundary ramified, the
ordered boundary pair remains intrinsic.  The `C=0` target is uniquely the
one carrying the generic irreducible residue-degree-two component (7.24)
inside the distinguished affine source.  Over the discriminant target, the
degree-two repeated-root cluster lies on the normalized boundary instead;
all generic affine source sheets are simple.  This affine-versus-boundary
distinction survives stabilization and replaces the ramification-only
ordering criterion used on the squarefree open.

Repeated roots satisfy `H(rho)=H'(rho)=0`, so their critical image can
coincide with the zero-cluster image.  More generally, several critical
parameters may have the same image, or another parameter may share the
root-one critical image.  These are conductor identifications on the target
critical curve.  They do not identify points on the normalized source
boundary edge.  The map-intrinsic package retains that normalization, while
the intersection with the unique root-one affine component still marks
`W=1`.  Consequently multi-double-root and affine-sheet image collisions do
not alter either marked principal ideal `(W)`, `(W-1)` or the divisor
recovered from (7.29).

> **Full collision theorem.**  Let `k` have positive characteristic and let
> `H` be an exact-degree normalized seed with an exact double zero:
>
> ```text
> H(0)=H'(0)=H(1)=0,       H'(1)=-1.
> ```
>
> Use whichever of the two complementary source charts is defined.  The
> complete intrinsic boundary package determines `H` exactly under stable
> polynomial left--right equivalence, with no squarefreeness, Hessian
> degree, critical-birationality, or collision-avoidance hypothesis.

The proof is now uniform.  The normalized critical source edge is an affine
`W`-line (possibly mapping non-birationally or inseparably to its target
image).  Its zero and root-one intersections fix `W` after stabilization.
Equation (7.29) recovers the effective divisor

```text
D_prim(H)=div(H/(W^2(W-1)))                           (7.31)
```

including residue extensions and multiplicities.  If `P_D` is the monic
polynomial with divisor `D_prim(H)`, then

```text
H=lambda W^2(W-1)P_D,
lambda=-1/P_D(1).                                    (7.32)
```

The denominator is nonzero because `1` is a simple distinguished root.
Equations (7.31)--(7.32) reconstruct `H` exactly.

This also answers the smaller-quotient question.  The full conductor and
self-intersection data are not required for seed recovery.  The stable
marked-edge quotient

```text
(A^1_W; (W), (W-1), D_prim(H))                       (7.33)
```

is map-intrinsic, survives every declared collision, and is already
coefficient-faithful.  Separate projectivized Hasse-channel divisors still
lose relative scales; (7.33) is the smaller exact divisor package.

## 8. Completed scope

The theorem completely removes the positive-characteristic integration
failure:

```text
p-typical Hasse coefficient data + H(0)
    <-> the bounded-degree polynomial H.              (8.1)
```

In particular, Frobenius directions no longer create
positive-dimensional fibers at the coefficient-data level.

Section 6 proves that the original derivative-only weighted map cannot have
an exact seed-recovery invariant.  Section 7 supplies the necessary
Frobenius enrichment and the complementary source chart, proves polynomial
Keller realization in every characteristic, and establishes stable
intrinsic recovery through clean, degree-drop, wild, repeated-root,
multi-double-root, and affine-sheet collision strata.  The exact marked
divisor quotient (7.33) is sufficient; the larger conductor package remains
useful for finer singularity questions but is not needed for seed recovery.

Within the declared normalized exact-double-seed problem, no
coefficient-recovery, realization, stable-intrinsic, or collision gap
remains.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_hasse_typical_seed_recovery.py
```

The checker verifies the reconstruction formula and exact kernel on monomial
bases for several primes and every bound through degree forty, verifies the
forced-channel minimality against every other Hasse order, and replays the
clean degree-eight `F_5` Hessian collision together with its separating
fifth Hasse channel.  It also verifies the five-member degree-twelve clean
family (6.6), equality of all ordinary-derivative construction data, direct
polynomiality and Jacobian one for their common weighted map, and separation
by `D^[5]`.  Finally it constructs all five enriched maps, checks their
polynomial correction, determinant, inverse-pencil identities, and the
marked transverse node at `c=2`.  The full equal-image Groebner bases certify
pairwise stable inequivalence of all five enriched maps.
The characteristic-two checks also verify the parity reconstruction and
singular old parameter, exhaust scalar reparameterizations through degree
eight as a regression of the general no-go proof, and check polynomiality,
the inverse identities, and Jacobian one for the replacement suspension
(7.15) in exact degrees three through twelve whenever the displayed
`F_2` normalization slice is nonempty.  Finally it verifies the radicial
factorization (7.22), compressed birational witnesses, and squarefree
primitive factors (7.26) in every degree from five through sixteen, the
symbolic quartic slice (7.27), and the absence of a non-birational compressed
map on every clean `F_2` seed through degree twelve.  The singular-parameter
quartic is also checked directly in characteristics three, five, and seven.
Finally, repeated-root collision seeds in characteristics two, three, and
five verify critical-image collision at `(0,0)` and exact reconstruction
from the marked primitive divisor with its full multiplicities.
