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

## 5. Determinant-one normalization

For a map with determinant `d`, translate its selected target to zero and
multiply one target coordinate by `d^(-1)`.  This is a polynomial target
automorphism, changes the determinant to one, preserves the complete fiber,
and does not change the stable left--right class.  All maps above may
therefore be presented with determinant one and the common fiber over the
origin.

## 6. Exact regression

Run

```bash
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
```

The checker verifies the all-degree tangent and boundary-clean identities,
the two common inverse pencils, the small rational quartic, the
quadratic-field cancellation specialization, the three common targets, and
the mod-`17` Rabin certificate.  The family-specific stable-boundary
calculations are independently checked by

```bash
.venv/bin/python scripts/verify_same_degree_stable_inequivalence.py
.venv/bin/python scripts/verify_quadratic_weighted_stable_separation.py
.venv/bin/python scripts/verify_quadratic_cancellation_intersection.py
```
