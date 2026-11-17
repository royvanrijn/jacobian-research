# Cross-resultants restrict every factorization, in all degrees

This note addresses factors without rational roots. It proves an
all-dimensional factor-count restriction, its plane quartic consequence,
and a precise form for the remaining two unequal factors. It does not
prove arbitrary-degree irreducibility or resolve either plane endpoint.
The registered claims are `HC4MYIF0/1/2` in
[MATH_STATUS.json](../../MATH_STATUS.json).

Let `n>=2`, let `K=(f_1,...,f_n):A^n->A^n` be complex Keller, and put

```
B=C[x_1,...,x_n], L=Frac(B), F=C(f_1,...,f_n),
P(T,z)=sum_(i=0)^m b_i(T)z^i, a=b_m!=0, m>=2,
S(z)=P(f_1,z)+f_2.
```

Every target coefficient, source degree and `m` is unrestricted.
The source field `L/F` is finite. We use the retained
[square-free pullback and divisor-image proof HC4MYHW0](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md).

## 1. Coprime leading factors give constant cross-resultants (HC4MYIF0)

The polynomial `S` is primitive over the UFD `B`. Indeed a source prime
dividing all its coefficients would divide a nonzero positive coefficient,
hence lie over `f_1=c`. Its divisibility of `b_0(f_1)+f_2` would contradict
dominance of that source divisor onto the target divisor `T=c`.

The hypersurface `S=0` in `A^(n+1)` is smooth. Its source partial
derivatives form the covector

```
basis coefficient row (P_T(f_1,z),1,0,...,0) times J(K).
```

It is nonzero at every point, because `J(K)` is invertible. Consequently
any distinct primitive polynomial factors of `S` have disjoint finite
zero loci: an intersection would make every partial derivative of their
product vanish. All factors are separable and have positive `z` degree.

Suppose a factorization over `L` is represented, by Gauss lemma, as

```
S=F_1(z)...F_k(z),  d_i=deg_z(F_i)>0,
l_i=leading_z(F_i),  product_i l_i=a(f_1).      (1)
```

Absorb the overall nonzero constant into one factor. Assume the `l_i`
are pairwise coprime in `B`. Then

```
Res_z(F_i,F_j)=lambda_ij in C*,  i!=j.         (2)
```

To see this, a source prime dividing this resultant cannot divide both
leading coefficients. Over the fraction field of the residue domain,
a zero resultant would therefore give a common finite root. Its closure
would intersect the two finite zero loci, impossible by smoothness.
The resultant has no source prime divisor and is thus a nonzero constant.

The leading-coefficient hypothesis matters. The two polynomials
`sz-1, sz-2` have disjoint finite zero loci, but resultant `-s`: their
leading coefficients share `s`. A finite-root argument alone misses
this common projective root at infinity. No Keller pair is asserted by
this control.

### 1.1 A coefficient test guaranteeing coprime leading factors

The condition

```
gcd(a,a',b_(m-1))=1                           (3)
```

guarantees pairwise coprime `l_i` in every source factorization (1).
Let `delta` divide two such leading coefficients, and let it lie over
`f_1=c`. Write `q=ord_c(a)`; square-free pullback gives
`ord_delta(a(f_1))=q`. If `q=1`, the product (1) cannot have two leading
coefficients divisible by `delta`. If `q>=2`, the coefficient of
`z^(m-1)` in (1) is a sum in which each term omits at most one `l_i`.
Every term is divisible by `delta`. Thus `b_(m-1)(c)=0`, contradicting
(3). This proof does not require a Newton polygon or splitting a
ramified cluster into separate roots.

In particular (3) covers every square-free `a` of arbitrary parameter
degree and every `a` coprime to `b_(m-1)`, with arbitrary multiplicities.

## 2. The valuation at V=infinity leaves at most two unequal factors

Under the coprime-leading hypothesis, no factorization (1) can have
`k>=3`, and a two-factor decomposition cannot have equal degrees.
This statement holds for every `n>=2` and `m>=2`.

Normalize by `U=a(T)z`:

```
Q(T,U,V)=a(T)^(m-1)(P(T,U/a(T))+V),
Q_i(U)=[a(f_1)^d_i/l_i] F_i(U/a(f_1)).
```

Each `Q_i` is monic and belongs to `B[U]`; all nonleading coefficients
have enough powers of `a(f_1)` to divide by `l_i`. Their product is
`Q(f_1,U,f_2)`. Resultant scaling and (2) give the exact identity

```
Res_U(Q_i,Q_j)
   =lambda_ij a(f_1)^(d_i*d_j)/(l_i^d_j l_j^d_i).       (4)
```

Take the discrete valuation of `F=C(f_1,f_3,...,f_n)(f_2)` at
`V=f_2=infinity`, normalized by `v(V)=-1`. Extend it to a common field
containing `L` and the roots of `Q`. Its restriction has `v(a(f_1))=0`.
Every difference of distinct roots of `Q` has valuation `-1/m`.

Here is a direct justification of the last statement. Set `V=tau^(-m)`
and `U=tau^(-1)X`. Then `tau^m Q` has reduction

```
X^m+a(T)^(m-1)
```

at `tau=0`, with `m` distinct nonzero roots over the algebraic closure
of the constant field. Each root lifts uniquely to a formal power
series: at each order its next coefficient is determined by a linear
equation with nonzero coefficient `m X(0)^(m-1)`. These are all `m`
roots. Differences have distinct nonzero leading coefficients, so their
valuations in the normalization `v(V)=-1` are exactly `-1/m`. No finite
degree search or general Puiseux theorem is used.

Put `ell_i=v(l_i)`. The left side of (4) is a product of `d_i*d_j`
cross-differences. Taking valuations gives, for every pair,

```
ell_i/d_i + ell_j/d_j = 1/m,
sum_i ell_i=0.                                (5)
```

If `k>=3`, any three equations make their three ratios equal to
`1/(2m)`. The remaining pair equations give the same value for every
other ratio. But then `sum ell_i=(sum d_i)/(2m)=1/2`, contradicting (5).
If `k=2` and `d_1=d_2=m/2`, (5) again gives `sum ell_i=1/2`, the same
contradiction.

Thus under (3), every reducible source polynomial has exactly two
irreducible factors of unequal positive degrees. More generally, without
(3), any three-factor or equal-two-factor decomposition requires a common
source prime among some leading coefficients; it lies over a common
zero of `a,a',b_(m-1)`. This is an obstruction to those factor patterns,
not to every unequal two-factor pattern.

## 3. Plane quartics close; arbitrary-degree survivors have two larger factors (HC4MYIF1)

Now take `n=2`. The retained
[HC4MYIE0 derivative-surface theorem](DERIVATIVE_POWER_MODIFICATION.md)
excludes rational roots under (3). Combining it with (5), a reducible
`P(f,Z)+g` under (3) must have exactly two irreducible factors of
degrees

```
2 <= d < e,     d+e=m.                        (6)
```

In particular every quartic satisfying

```
gcd(a,a',b_3)=1
```

stays irreducible under every plane Keller base change. A proper quartic
factorization either has a linear factor or splits `2+2`; both are
excluded. This closes all square-free-leading quartics of arbitrary
parameter degree, and all quartics with `gcd(a,b_3)=1`, with arbitrary
leading multiplicities. No nilpotent root-group, critical-discriminant,
or special Keller-component fibre hypothesis is imposed. The conclusion
is full quartic irreducibility, not just absence of rational roots.

Consequently every closed stable embedding `A^(r+2)->A^(r+3)`, `r>=1`,
carrying a plane Keller pair through `(T,xy-P)` in this quartic class
forces the pair invertible. The retained
[HC4MYHZ2 filtered output proof](LEADING_FIBRE_UNIT_GATE.md) applies:
it supplies the filtered basis, all affine height-one geometric fibre
conditions, the generic root lattice and the global component pencil.
Root transitivity suffices; full splitting-group preservation is not
asserted.

For higher degrees, for example, the `2+3` pattern in degree five and
the `2+4` pattern in degree six are not excluded by (5). The equal and
three-factor exclusions are all-degree assertions; the quartic corollary
does not imply full irreducibility for `m>=5`.

## 4. The unequal case has an explicit resultant power and a nonnormal order (HC4MYIF2)

For any `n>=2`, suppose (1) has two factors of degrees `d<e` and coprime
leading coefficients. Write their leading coefficients as `p,D`, so
`a(f_1)=pD`, and put `r=e-d`, `alpha=e(d-1)`. With `Q_d,Q_e` monic as
above, (4) becomes

```
R=Res_U(Q_d,Q_e)=lambda a(f_1)^alpha D^r.       (7)
```

The factors' coefficients and `R` are integral over the target ring
`A=C[f_1,...,f_n]`: they are symmetric polynomials in roots of the monic
target equation `Q`. In the linear-factor case `d=1`, `R` is the
normalized derivative at that root, and (7) is exactly its
`(m-2)`nd-power form, up to a nonzero constant. Equation (7) is its
all-factor replacement.

Solving (5) gives, at every extension of the target infinity valuation,

```
v(p)=d e/[m(e-d)],    v(D)=-d e/[m(e-d)].       (8)
```

The leading coefficient of the larger factor has a pole at every such
place, whereas that of the smaller factor has a zero. These rational
values solve the two-factor equations; they do not themselves give a
contradiction.

They do give an exact field-degree constraint:

```
m(e-d)/gcd(d e,m(e-d)) divides [L:F].          (9)
```

Indeed every `F`-conjugate of `D` still satisfies the conjugate of (7),
whose cross-resultant has valuation `-de/m`. Thus every conjugate has
the same valuation (8). The valuation of the norm from `L` to `F` is
`-[L:F]de/[m(e-d)]`, which must be an integer. For `d=2,e=3`, (9)
requires a multiple of five. This is a conditional restriction on all
Keller degrees, not a search over particular source degrees.

There is also an explicit nonnormal integral order. Put

```
N=ceil(alpha/r),     zeta=a(f_1)^N D.
```

Then `zeta^r=lambda^(-1)a(f_1)^(Nr-alpha) R` is integral over `A` by
(7), so `zeta` is integral over `A` by transitivity of integrality.
It is a polynomial on the source, with its zero divisors supported only
over `a=0`.

Moreover `zeta` is not in `F`. To justify this last point, the retained
square-free/coprime divisor pullbacks imply `B intersect F=A`: a target
rational function with a denominator prime cannot pull back to a source
polynomial, since that prime's nonconstant pullback is coprime to the
numerator pullback. If `D` belonged to `F`, both `D` and `p=a/D` would
therefore belong to `A`. Their product `a(f_1)` forces both into
`C[f_1]`, contradicting the nonzero valuations (8). Multiplication by a
power of `a(f_1)` does not change membership in `F`.

Consequently the finite monogenic order `A[zeta]` is **nonnormal**.
If it were normal, the retained
[HC4MYHY0 normal monogenic theorem](NORMAL_MONOGENIC_COVERS.md) would
preserve irreducibility of its monic minimal polynomial under Keller
base change; that polynomial already has source root `zeta` and has
degree greater than one. This contradiction proves nonnormality.
The finite normalization exists and embeds in `B`, but it is not thereby
proved globally monogenic. No purity argument is applied to this
nonnormal primitive model.

## 5. Boundary of the result, literature audit and replay

The new factor-count obstruction uses only primitive factorization,
smoothness, resultant scaling, the infinity valuation and constant source
units. It does not import a new classification theorem. Section 4 reuses
the explicitly sourced normal-monogenic theorem for its nonnormality
conclusion.

The alternative ruled-surface route was inspected in primary literature.
Miyanishi's 2007 Theorems 1.1/1.2 in *Q-factorial subalgebras of a
polynomial ring* impose Q-factoriality, and the 2005 affine pseudo-covering
discussion imposes almost surjectivity when invoking that notion.
Neither property follows from a dominant map to the normalized
`P(T,z)+V` surface. Neither classification is imported. The
[source audit](sources/factor-resultant-context.json) records the inspected
statements and the missing hypotheses.

Under (3), the higher-degree residual is now precisely (6), with (7)--(9)
and an explicit nonnormal integral order. When (3) fails, common source
leading factors and actual intermediate pole behavior remain separate
obstructions. Other output geometries, singleton-zero and unit-ideal
plane endpoints and unrestricted HC4/JC2 remain open. None of the stated
factor patterns is asserted realizable by a Keller pair.

Run [the exact replay](verify_factor_resultants.py):

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_factor_resultants.py
```

The [receipt](receipt-factor-resultants.json) contains finite identity and
boundary controls. They accompany the unrestricted proofs above. The
[validation](validation-factor-resultants.json) records scoped checks
separately from whole-repository failures. The automorphic actual
unit-Hessian embedding controls from the retained
[derivative note](DERIVATIVE_POWER_MODIFICATION.md) remain valid for every
`P`; no nonlinear Hessian invariance is used. Formal verification,
independent replay, external review and novelty priority are not claimed.
