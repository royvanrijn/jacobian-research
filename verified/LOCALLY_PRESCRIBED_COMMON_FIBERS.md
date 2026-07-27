# Locally prescribed fibers in one fixed inequivalent pair

The [common-arithmetic-fiber theorem](COMMON_ARITHMETIC_FIBERS.md) gives, for
every `N>=4`, two fixed stably inequivalent Keller maps sharing the pencil

\[
 P_{N,u}(T)=T^N+T^3-2T^2+T+u.                         \tag{1}
\]

The [local-to-global theorem](LOCAL_GLOBAL_KELLER_FIBERS.md) prescribes
arbitrary finite local data by allowing the map to depend on the resulting
polynomial.  This note keeps both maps fixed.  The necessary price is a
family-relative compatibility condition: every requested local algebra must
occur at some local parameter in the one-dimensional pencil (1).

## 1. The fixed determinant-one pair

Put

\[
 H_N(T)=T^N+T^3-2T^2,\qquad G_N(T)=H_N(T)+T.
\]

The weighted map `F_N^wt` attached to `H_N` has determinant `1-N`; at

\[
                         q_u^{\rm wt}=
                         \left({u\over1-N},-1,1\right)             \tag{2}
\]

its inverse polynomial is (1).  The quadratic-gauge map `F_N^quad` attached
to `G_N` has determinant `-2`; at

\[
                         q_u^{\rm quad}=(1,0,-2u)                 \tag{3}
\]

its inverse polynomial is again (1).  These maps are fixed as `u` varies and
are stably polynomially left--right inequivalent.

Apply the fixed target scalings

\[
 (A,B,C)\longmapsto\left({A\over1-N},B,C\right),\qquad
 (\Pi,B,C)\longmapsto(\Pi,-B/2,C).                                \tag{4}
\]

The resulting two fixed maps have determinant one.  Their common targets
are

\[
 \widetilde q_u^{\rm wt}=
 \left({u\over(1-N)^2},-1,1\right),\qquad
 \widetilde q_u^{\rm quad}=(1,0,-2u).                             \tag{5}
\]

For every squarefree `P_(N,u)`, both fibers are complete and isomorphic to
`Spec Q[T]/(P_(N,u))`.

## 2. Family-relative local specifications

Fix a finite set `S` of rational primes.  At each `p in S`, prescribe a
rank-`N` finite étale algebra `A_p/Q_p` together with a local witness
`u_p in Q_p` such that

\[
                         A_p\simeq
             \mathbb Q_p[T]/(P_{N,u_p}).                          \tag{6}
\]

This witness condition is load-bearing.  A one-parameter pencil cannot
realize every rank-`N` local algebra.

Choose also a nonempty real interval `I_infty` on which `P_(N,u)` is
squarefree with a fixed number `r_1` of real roots.  Further unramified
Frobenius partitions may be included by choosing finite-field parameter
witnesses and regarding their residue classes as additional local
conditions.

> **Locally prescribed fixed-common-fiber theorem.**  Under these
> family-relative compatibility conditions, there are infinitely many
> `u in Q intersect I_infty` such that
>
> \[
> \mathbb Q[T]/(P_{N,u})\otimes\mathbb Q_p\simeq A_p
> \quad(p\in S),                                                    \tag{7}
> \]
>
> `P_(N,u)` has `r_1` real roots, and `P_(N,u)` is irreducible over `Q`.
> Consequently, the two fixed determinant-one, stably inequivalent maps in
> Section 1 share infinitely many connected complete fibers with all the
> prescribed local data.

### Proof

For each `p`, the isomorphism class in (6) is open in the coefficient
topology.  Because only the constant coefficient varies, its inverse image
is an open neighborhood `Omega_p` of `u_p` in `Q_p`.  Weak approximation
meets

\[
                         I_\infty\times\prod_{p\in S}\Omega_p.     \tag{8}
\]

The two common-fiber identities are polynomial identities in `u`, so every
rational point of (8) gives the same prescribed complete fiber in both
fixed maps.  Finally, `P_(N,U)` is irreducible over `Q(U)` because it is
primitive and linear in `U`.  Hilbert irreducibility with weak approximation
gives infinitely many irreducible rational specializations inside (8).
Their real root count is constant on `I_infty`.  This proves the theorem.
\(\square\)

### Constructive connectedness

Suppose, in addition, that there are a prime `ell` and a residue
`bar u_ell` for which `P_(N,bar u_ell)` is irreducible over `F_ell`.  Add

\[
                              u\equiv\bar u_\ell\pmod\ell          \tag{9}
\]

to (8).  Every rational lift with `ell`-integral coefficients is then
irreducible over `Q`.  A denominator-compatible CRT grid gives infinitely
many such lifts in `I_infty`, so this variant is fully constructive and does
not use Hilbert irreducibility.

## 3. Automatic parameter radii

Assume that a local witness `u_p` is rational and `p`-integral.  Put

\[
 D_p=v_p(\operatorname{Disc}P_{N,u_p}),\qquad
 m_p=2D_p+1.                                                       \tag{10}
\]

The automatic local-stability theorem gives

\[
 v_p(u-u_p)\ge m_p
 \quad\Longrightarrow\quad
 \mathbb Q_p[T]/(P_{N,u})\simeq
 \mathbb Q_p[T]/(P_{N,u_p}),                                     \tag{11}
\]

because `P_(N,u)-P_(N,u_p)=u-u_p`.  Thus full local-algebra preservation
reduces to one scalar congruence.  The function

```python
jcsearch.common_fibers.synthesize_fixed_common_fiber_parameter
```

derives all exponents (10), performs the prime-power CRT, meets a rational
real interval, and returns the parameter, polynomial, and local
certificates.

The discriminant exponent is conservative.  A factorwise analysis can
replace `m_p` by a smaller certified value without changing the theorem.

## 4. A ramified sextic common fiber

Take `N=6`.  Prescribe the local centers

\[
                         u_2=-1,\qquad u_3=2,\qquad u_5=1.         \tag{12}
\]

At `2`,

\[
 P_{6,-1}=(T-1)(T^2+1)(T^3+T^2+1).                              \tag{13}
\]

The quadratic defines the ramified extension `Q_2(i)`, while the cubic is
irreducible modulo `2` and defines the unramified cubic `U_(2,3)`.  Hence

\[
 \mathbb Q_2[T]/(P_{6,-1})
 \simeq\mathbb Q_2\times\mathbb Q_2(i)\times U_{2,3}.             \tag{14}
\]

Its polynomial discriminant is `4464`, of `2`-adic valuation four, so the
automatic parameter precision is nine.

At `3`, reduction gives

\[
 P_{6,2}\equiv
 (T+2)^2(T^2+1)(T^2+2T+2)\pmod3.                                 \tag{15}
\]

The two quadratics are irreducible and lift to two copies of the unramified
quadratic extension.  After `T=1+X`, the coefficient valuations of
`P_(6,2)(1+X)` in ascending order are

\[
                         (1,1,0,1,1,1,0).                         \tag{16}
\]

The first Newton segment has length two and slope `-1/2`, giving one
ramified quadratic factor.  Thus the full local algebra is one ramified
quadratic times two unramified quadratics.  Its discriminant valuation is
one, so the automatic precision is three.

Modulo `5`, the polynomial `P_(6,1)` is irreducible.  Its discriminant is a
`5`-adic unit, so precision one both preserves the unramified degree-six
field and certifies global connectedness.

Use the real interval

\[
                              {1\over2}<u<{3\over2}.               \tag{17}
\]

The parameter CRT modulus and its first denominator-one lift are

\[
 2^9\,3^3\,5=69120,\qquad
 \boxed{u={95231\over69121}}.                                    \tag{18}
\]

The resulting polynomial

\[
 \boxed{
 P(T)=T^6+T^3-2T^2+T+{95231\over69121}
 }                                                               \tag{19}
\]

is irreducible modulo `5`, has exactly two real roots, and therefore defines
a sextic field of signature `(2,2)`.  Equations (10)--(11) identify its
completions at `2` and `3` with the ramified algebras above.

The two fixed determinant-one maps have common targets

\[
 \boxed{
 \widetilde q^{\rm wt}
   =\left({95231\over1728025},-1,1\right),\qquad
 \widetilde q^{\rm quad}
   =\left(1,0,-{190462\over69121}\right).
 }                                                               \tag{20}
\]

They are stably inequivalent, while both complete fibers are
`Spec Q[T]/(P)`.

## Verification

Run

```bash
.venv/bin/python scripts/verify_locally_prescribed_common_fibers.py
```

The checker verifies the local factorizations, Newton polygon, discriminant
radii, parameter CRT, exact real-root count, inert-prime irreducibility,
both fixed inverse equations, determinant-one target transport, and the two
common targets.
