# Common arithmetic fibers of stably inequivalent Keller maps

## 1. The invariant

Let `K` be a characteristic-zero field and let `A/K` be a finite etale
algebra.  Write

\[
 \mathcal R_K(A)=
 \left\{
 \begin{array}{c}
 \text{stable polynomial left--right classes of Keller maps}\\
 \text{having a complete fiber isomorphic to }\operatorname{Spec}A
 \end{array}
 \right\}.
\]

Left--right equivalence and stabilization are taken over `K`.  A complete
fiber means a regular fiber whose length equals the geometric degree of the
map.  This convention makes the property invariant under left--right
equivalence and under adjoining identity variables.

## 2. Two fixed maps over `Q`

For every `N>=4`, put

\[
 H_N(T)=T^N+T^3-2T^2,\qquad
 P_{N,u}(T)=H_N(T)+T+u.
\]

The polynomial `P_(N,u)` is tangent-admissible for every `u`:

\[
 P_{N,u}(1)-P_{N,u}(0)=P_{N,u}'(0)=1,
\]

\[
 P_{N,u}'(0)-P_{N,u}'(1)=1-N\ne0,
\]

\[
 P_{N,u}''(1)=N(N-1)+2\ne2(N-1).
\]

Moreover,

\[
 H_N(T)=T^2(T^{N-2}+T-2)
\]

is boundary-clean.  Indeed zero is an exact double root.  If a nonzero
common root of `T^(N-2)+T-2` and its derivative existed, the derivative
equation and the original equation would force

\[
 T=\frac{2(N-2)}{N-3}>0,\qquad
 T^{N-3}=-\frac1{N-2},
\]

which is impossible.

The weighted map attached to `H_N` is fixed as `u` varies.  It has
determinant `1-N`, and at

\[
 q_u^{\rm wt}=\left(\frac{u}{1-N},-1,1\right)
\]

its inverse polynomial is `P_(N,u)`.

The root-engineered quadratic-gauge map attached to

\[
 G_N(T)=H_N(T)+T
\]

is also fixed as `u` varies.  Its determinant is `-2`, and at

\[
 q_u^{\rm quad}=(1,0,-2u)
\]

its inverse polynomial is again `P_(N,u)`.

Whenever `P_(N,u)` is squarefree, both complete fibers are isomorphic to

\[
 \operatorname{Spec}\mathbb Q[T]/(P_{N,u}).
\]

The polynomial `P_(N,U)` is irreducible over `Q(U)`: it is primitive and
linear in `U`, so any factorization in `Q[T,U]` would have a factor
independent of `U`, and comparison of the coefficient of `U` makes that
factor a unit.  Hilbert irreducibility therefore supplies a Hilbert subset
of rational `u` for which `P_(N,u)` is irreducible.

The boundary-clean weighted map and every admissible quadratic-gauge map of
the same degree are stably inequivalent: after deleting the intrinsic second
boundary vertex, the normalized ramified target strata are

\[
 \mathbb A^1\times\mathbb G_m
 \quad\text{and}\quad
 \mathbb G_m^2,
\]

with unit ranks one and two.  Thus:

> **Fixed-pair common-fiber theorem.**  For every `N>=4`, two fixed
> `Q`-defined Keller maps of geometric degree `N` share complete fibers
> \[
> \operatorname{Spec}\mathbb Q[T]/(P_{N,u})
> \]
> for every squarefree specialization `u`.  The two maps are not stably
> polynomially left--right equivalent.  For a Hilbert subset of rational
> `u`, the common fiber is connected.  In particular,
> \[
> |\mathcal R_{\mathbb Q}(\mathbb Q[T]/(P_{N,u}))|\ge2.
> \]

The family-relative adelic strengthening is proved in
[locally prescribed common fibers](LOCALLY_PRESCRIBED_COMMON_FIBERS.md):
whenever selected local algebras occur at local parameters of this pencil,
infinitely many connected common fibers realize them simultaneously.
Automatic discriminant radii make the parameter congruences explicit.

## 3. The small rational quartic

The particularly small quartic

\[
 P(T)=2T^4-T^3-T^2+T+1
\]

also lies in the overlap.  Its weighted seed is

\[
 H(T)=P(T)-1-T=T^2(T-1)(2T+1).
\]

The weighted determinant and target are

\[
 -3,\qquad \left(-\frac13,-1,1\right).
\]

The quadratic-gauge seed is

\[
 G(T)=P(T)-1,
\]

and its target is `(1,0,-2)`.  Reduction modulo `3` proves that `P` is
irreducible, and its discriminant is `1556`.  Hence

\[
 A=\mathbb Q[T]/(P)
\]

is a quartic field with `|R_Q(A)|>=2`.

## 4. Three fixed maps over `Q(sqrt(-2))`

Let

\[
 K=\mathbb Q(\eta),\qquad \eta^2=-2,
\]

and define the fixed quartic pencil

\[
 \Psi_R(T)=
 T-\frac92T^2+(8+2\eta)T^3
 -\left(\frac72+2\eta\right)T^4-R.
\]

It is the type-`(2,1)` cancellation inverse at target

\[
 (\Pi,Q,R)=(4+\eta,3,R).
\]

The associated polynomial cancellation map is fixed: take

\[
 \theta=2+\eta,\qquad
 h(A)=\theta+(4\theta-6)A.
\]

The relation `theta^2-4theta+6=0` is exactly its polynomiality condition,
and the map has determinant `-1`.

The same pencil is tangent-admissible:

\[
 \Psi_R(1)-\Psi_R(0)=\Psi_R'(0)=1,
\]

\[
 c=\Psi_R'(0)-\Psi_R'(1)=-1+2\eta,
\]

\[
 \Psi_R''(1)-2(\Psi_R'(1)-\Psi_R'(0))=-5-8\eta.
\]

Its fixed weighted seed is

\[
 H(T)=\Psi_R(T)-\Psi_R(0)-T
 =T^2(T-1)
 \left(
 \frac92-\left(\frac72+2\eta\right)T
 \right),
\]

which is boundary-clean.  The weighted target is

\[
 q_R^{\rm wt}=\left(\frac{-R}{-1+2\eta},-1,1\right).
\]

The fixed quadratic-gauge seed is

\[
 G(T)=\Psi_R(T)-\Psi_R(0)=\Psi_R(T)+R,
\]

with target

\[
 q_R^{\rm quad}=(1,0,2R).
\]

The polynomial `Psi_R` is irreducible over `K(R)` by the same
degree-one-in-`R` argument.  Hilbert irreducibility supplies infinitely many
`R in K` for which it is irreducible.

The three fixed maps are pairwise stably inequivalent.  Unit rank separates
the weighted map from the two reciprocal constructions.  On the common
ramified torus, the cancellation relative Fitting generator has affine
Laurent-support rank one, while the quadratic-gauge generator has rank two.
Independently, the intrinsic boundary-contact nilpotency indices are

\[
 1,\qquad6,\qquad2
\]

for weighted, cancellation, and quadratic gauge.

> **Fixed-triple common-fiber theorem.**  Three fixed Keller maps over
> `K=Q(sqrt(-2))`, all of geometric degree four, share the complete fiber
> \[
> \operatorname{Spec}K[T]/(\Psi_R)
> \]
> for every squarefree specialization `R`.  They are pairwise stably
> polynomially left--right inequivalent.  On a Hilbert subset of `K`, the
> common fiber is connected, so
> \[
> |\mathcal R_K(K[T]/(\Psi_R))|\ge3.
> \]

At `R=-1`, irreducibility has a short finite-field certificate.  Reduce
`2 Psi_(-1)` at `(17,eta-10)`.  After monic normalization the result is

\[
 f(T)=T^4+14T^3+2T^2+9T+9\in\mathbb F_{17}[T].
\]

One has

\[
 T^{17^4}-T\equiv0\pmod f,
\]

\[
 T^{17^2}-T\equiv6T^3-T^2-6T+3\pmod f,
\qquad
 \gcd(f,T^{17^2}-T)=1.
\]

Rabin's criterion proves that `f`, and hence `Psi_(-1)`, is irreducible.

## 5. One small rational-coefficient field in all three maps

The Hilbert-family statement has the following small individual witness.
Keep

\[
 K=\mathbb Q(\eta),\qquad \eta^2=-2,
\]

and put

\[
 \boxed{p(W)=9W^4-19W^3+10W^2-8W-4.}                 \tag{5.1}
\]

Although the cancellation map requires `K`, the polynomial (5.1) has
rational integer coefficients.  In fact it defines a quartic field `L/Q`
and remains irreducible over `K`, so

\[
 \boxed{A=K[W]/(p)=K\mathbin{\mathop{\otimes}_{\mathbb Q}}L}
                                                               \tag{5.2}
\]

is a connected quartic finite etale `K`-algebra.

Here are the three fixed maps and their selected targets.

### 5.1 Weighted

Take the integral boundary-clean seed

\[
 H(W)=9W^4-19W^3+10W^2=W^2(W-1)(9W-10).
\]

It has `c=-H'(1)=1`, `H''(1)/c=14`, and weighted parameter
`a_0=-15/16`.  To display the actual determinant-one polynomial map, put

\[
 u=1+xy,\qquad
 \gamma=1-\frac{15}{16}xy+x^2z,\qquad W=u\gamma,
\]

\[
 p_H(W)=W(36W^2-57W+20),\qquad
 q_H(W)=W^2(27W^2-38W+10),
\]

and set

\[
 F^{\rm wt}=
 \left(
  \frac{u+q_H(W)/\gamma^2}{x^2},\
  \frac{1+p_H(W)/\gamma}{x},\
  x\gamma
 \right).                                             \tag{5.3}
\]

The two displayed quotients are polynomials by weighted admissibility, and
`\det DF^{wt}=1`.  At the target

\[
 \boxed{q^{\rm wt}=(-4,8,1)}
\]

the inverse equation is exactly

\[
 H(W)-8W-4=p(W).                                      \tag{5.4}
\]

### 5.2 Cancellation

Take the fixed type-`(2,1)` map associated to

\[
 \theta=2+\eta,\qquad h(A)=\theta+(4\theta-6)A.
\]

Its complete denominator-free formula is
[displayed here](SAME_DEGREE_STABLE_INEQUIVALENCE.md#12-a-cancellation-map).
Write its target coordinates as `(\Pi,Q,R)`.  At

\[
 \boxed{
 q^{\rm can}
 =\left(\frac4{11},1,-\frac{22481}{23232}\right)
 }
                                                               \tag{5.5}
\]

make the affine generator change

\[
 T=\frac14+3W.
\]

The cancellation inverse polynomial then satisfies the exact identity

\[
 \boxed{
 \Psi_{q^{\rm can}}\left(\frac14+3W\right)
 =-\frac{36}{121}p(W).
 }                                                        \tag{5.6}
\]

The map has determinant `-1`.  Differentiating (5.6) shows that every
simple root of `p` has nonzero cancellation reconstruction denominator.

### 5.3 Quadratic gauge

Use the integral seed

\[
 G(S)=9S^4-19S^3-8S.
\]

The quadratic coefficient of `p` has been sheared into the selected target;
this removes the `tq` term from the map.  The explicit root-engineered map
is obtained by putting

\[
 t=1+xy,\qquad
 q=t^2z+\frac8{19}y^2(1+3t)
\]

and

\[
 \begin{aligned}
 F^{\rm quad}=\biggl(&tq,\\
 &y+\frac{57}{8}xq-\frac92t^2x^2q^4,\\
 &x(5-3t)-\frac{19}{8}x^3z
       +\frac94(xq)^4\biggr).                          \tag{5.7}
 \end{aligned}
\]

It has determinant `-2`.  At

\[
 \boxed{q^{\rm quad}=\left(1,\frac52,-1\right)}
\]

its inverse equation is

\[
 G(S)-\frac{-8}{2}\left(\frac52S^2-1\right)=p(S).
                                                               \tag{5.8}
\]

The expanded support counts of `F^wt` and `F^quad` are respectively
`(16,14,3)` and `(7,51,38)`, with component degrees `(12,11,4)` and
`(7,26,24)`.  The quadratic shear used above removes seven terms from its
second component.  The cancellation map is fixed independently of the
selected target.

### 5.4 Connectedness, completeness, and stable fingerprints

Modulo `17`, monic normalization of (5.1) is

\[
 f(W)=W^4-4W^3+3W^2+W-8.
\]

One has `7^2=-2 mod 17`, so the chosen prime splits in `K`.  Moreover,

\[
 W^{17^4}-W\equiv0\pmod f,
\]

\[
 W^{17^2}-W
 \equiv2W^3+8W^2+W-8\pmod f,\qquad
 \gcd(f,W^{17^2}-W)=1.
\]

Rabin's criterion proves irreducibility over the residue field and hence
over `K`.  Thus `p` is squarefree.  Equations (5.4), (5.6), and (5.8),
together with the three reconstruction theorems, show that all three
selected fibers are complete and isomorphic to `Spec A`.

The three one-line stable fingerprints are

\[
 \boxed{\text{unit ranks }(1,2,2)}
 \quad
 \text{for }({\rm wt},{\rm can},{\rm quad}),             \tag{5.9}
\]

\[
 \boxed{\text{reciprocal Fitting Laurent ranks }(1,2)}
 \quad
 \text{for }({\rm can},{\rm quad}),                     \tag{5.10}
\]

\[
 \boxed{\text{boundary nilpotency indices }(1,6,2)}
 \quad
 \text{for }({\rm wt},{\rm can},{\rm quad}).             \tag{5.11}
\]

Unit rank separates the weighted map from both reciprocal maps;
Laurent-support rank separates cancellation from quadratic gauge; and
(5.11) independently separates all three.

The affine-generator search leading to (5.1) is exact and bounded.  For the
quartic cancellation inverse, the weighted tangent-chord constraint reduces
to

\[
 y^2=-18a^2+24a-2.
\]

Enumerating the rational parameterization through `(a,y)=(1,2)`, with
`|num(k)|<=4`, `den(k)<=4`, `|num(rho)|<=12`, `den(rho)<=6`, and constant
term of the linear-normalized polynomial in `{1,-1,1/2,-1/2}`, selects
`a=1/11`, chord step `12/11`, and scale `11/4`.  These give
`T=1/4+3W`.  Among all certified candidates in that box, (5.1) uniquely
minimizes first the largest and then the sum of the absolute primitive
coefficients.

This is a bounded presentation-minimality certificate, not a global
minimality theorem under arbitrary affine generators or polynomial
left--right changes.

## 6. Determinant-one normalization

For a map with determinant `d`, translate its selected target to zero and
multiply one target coordinate by `d^(-1)`.  This is a polynomial target
automorphism, changes the determinant to one, preserves the complete fiber,
and does not change the stable left--right class.  All maps above may
therefore be presented with determinant one and the common fiber over the
origin.

## 7. Exact regression

Run

```bash
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
.venv/bin/python scripts/search_cross_family_collision.py
```

The checker verifies the all-degree tangent and boundary-clean identities,
the two common inverse pencils, the small rational quartic, the
quadratic-field cancellation specialization, the three common targets, and
both mod-`17` Rabin certificates.  The second command reproduces the bounded
affine-generator search and its coefficient-minimal output.  The
family-specific stable-boundary calculations are independently checked by

```bash
.venv/bin/python scripts/verify_same_degree_stable_inequivalence.py
.venv/bin/python scripts/verify_quadratic_weighted_stable_separation.py
.venv/bin/python scripts/verify_quadratic_cancellation_intersection.py
```
