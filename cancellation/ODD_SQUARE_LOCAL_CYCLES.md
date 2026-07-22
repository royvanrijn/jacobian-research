# Odd square family: isolated tame cycles

This note gives a structural Galois upgrade for the square-discriminant
family

\[
 m=2a^2-1,\qquad r=1,\qquad a\text{ odd}.
\]

It does not yet prove that infinitely many members have Galois group
`A_m`.  It does prove `A_m`, conditional only on irreducibility, on an
explicit arithmetic subfamily, and it identifies the two separate open
inputs needed to make that subfamily infinite.

Put

\[
 L=m+2=2a^2+1
\]

and use the reciprocal geometric-derivative polynomial

\[
 f_m(x)=x^m+2x^{m-1}+\cdots +(m+1).
\]

This is [Ishida's polynomial](https://arxiv.org/abs/1701.01160)
`f_(1,m+1)`.  It has the same splitting field
and the same natural Galois action as the cancellation polynomial
`M_(m,1)`.  The key identity is

\[
 \boxed{(x-1)^2f_m(x)=x^L-Lx+(L-1).}                 \tag{1}
\]

The discriminant is a rational square for odd `a`, so an irreducible member
has transitive Galois group `G<=A_m`.  The remaining task is to prove
`A_m<=G`.

## 1. A prime factor of `L` isolates a cycle

### Theorem 1 (isolated tame cycle)

Suppose

\[
 L=sq,\qquad q\ge5\text{ prime},\qquad q\nmid s.
\]

Then the Galois group of `f_m` contains a pure `(q-2)`-cycle: one cycle of
length `q-2`, fixing all other roots.

### Proof

Write the right side of (1) as

\[
 T(x)=x^L-Lx+(L-1).
\]

Work over an unramified extension of `Q_q` containing the `s`-th roots of
unity.  Modulo `q`,

\[
 T(x)\equiv x^{sq}-1=(x^s-1)^q.                     \tag{2}
\]

At the residue root `1`, put `x=1+y` and divide by the known factor `y^2`.
The coefficient of `y^i` in `f_m(1+y)` is
`binom(L,i+2)`.  Since `q || L`,

\[
 v_q\binom Lj=1\quad(2\le j<q),
 \qquad
 v_q\binom Lq=0.                                    \tag{3}
\]

Thus the Newton polygon begins with the single edge

\[
 (0,1)\longrightarrow(q-2,0).                       \tag{4}
\]

It gives a tamely and totally ramified root cluster of size `q-2`.

For every other `s`-th root of unity `zeta`, expand `T(zeta+y)`.  Its
constant and linear coefficients have valuation one,

\[
 T(\zeta)=L(1-\zeta),\qquad
 T'(\zeta)=L(\zeta^{-1}-1),
\]

while the same binomial calculation gives the edge

\[
 (0,1)\longrightarrow(q,0).                         \tag{5}
\]

Hence every other ramified cluster has size `q`.

It remains to check support, rather than merely the order of inertia.  Let
`P` be wild inertia in the local splitting field.  It fixes the
`(q-2)`-cluster pointwise, and tame inertia acts on that cluster as a
`(q-2)`-cycle.  On a `q`-cluster, normality forces the `P`-orbits to be
either singletons or the whole cluster.  The singleton case would make a
prime-to-`q` cyclic quotient transitive on `q` points, which is impossible.
Thus `P` acts regularly as `C_q`.  The inertia image on this cluster lies in
the normalizer `AGL_1(q)` of that regular group, so the image of a tame
complement has order dividing `q-1`.

Choose a generator `tau` of a tame inertia complement.  Its image on the
first cluster generates `C_(q-2)`.  The element `tau^(q-1)` still generates
that image because `gcd(q-2,q-1)=1`, and it fixes every `q`-cluster by the
preceding order bound.  Therefore `tau^(q-1)` is a pure `(q-2)`-cycle.  QED

The last paragraph is the support argument missing from the general
fixed-row Newton construction: here the exceptional cluster has tame size
`q-2`, while every competing cluster has prime size `q` and tame normalizer
of order dividing `q-1`.

## 2. The exact primitivity test

### Lemma 2 (one-cycle block test)

Let a transitive group of degree `m` contain a pure `k`-cycle.  If it
preserves blocks of size `c`, with `1<c<m`, then

\[
 c\ge k\qquad\text{or}\qquad c\mid k.                \tag{6}
\]

If the cycle fixes every block setwise, its support lies in one block and
`c>=k`.  If it moves `ell>1` blocks, no fixed point can lie in a moved block.
The union of those blocks is therefore exactly the support of the one
nontrivial cycle, so `k=ell*c` and `c|k`.

Let `lambda(m)` be the least prime divisor of the odd integer `m`.  Its
largest proper divisor is `m/lambda(m)`.  Theorem 1 and Lemma 2 therefore
give the following criterion.

### Corollary 3 (conditional `A_m` criterion)

Assume that `f_m` is irreducible.  If some prime `q>=5` satisfies

\[
 q\parallel L,\qquad 1< s=L/q,
\]

and

\[
 \gcd(m,q-2)=1,
 \qquad
 q-2>\frac{m}{\lambda(m)},                            \tag{7}
\]

then

\[
 \operatorname{Gal}(f_m/\mathbb Q)=A_m.              \tag{8}
\]

Indeed, (7) makes both alternatives in (6) impossible, so transitivity
upgrades to primitivity.  The `(q-2)`-cycle fixes `(s-1)q>=3` roots.  Jones's
classification theorem for primitive groups containing a cycle then gives
`A_m<=G`; see
[Jones, Corollary 1.3](https://arxiv.org/abs/1209.5169).  The square
discriminant gives `G<=A_m`.

The coprimality in (7) can be written without `m`:

\[
 \gcd(m,q-2)=\gcd(s-1,q-2),                           \tag{9}
\]

because `m=sq-2`, `q-2` is odd, and
`sq-2=2(s-1) mod (q-2)`.

## 3. A clean prime-cofactor subfamily

If

\[
 2a^2+1=3q\qquad(q\ge5\text{ prime}),                \tag{10}
\]

then (7) is automatic.  Here `m=3q-2`,
`gcd(m,q-2)=gcd(2,q-2)=1`, and

\[
 \frac{m}{\lambda(m)}\le\frac m3<q-2.
\]

Consequently:

> If `a` is odd, `(2a^2+1)/3` is prime, and `f_(2a^2-1)` is irreducible,
> then its Galois group is `A_(2a^2-1)`.

Writing `a=6t+1` or `6t-1`, the prime in (10) is respectively

\[
 q=24t^2+8t+1,
 \qquad
 q=24t^2-8t+1.                                       \tag{11}
\]

These primitive quadratics have no fixed prime divisor, so the usual
Bunyakovsky heuristic predicts infinitely many prime values.  Infinitude is
not known.  Moreover, irreducibility of `f_m` is itself a special case of
the geometric-derivative irreducibility conjecture.  Thus (10) is a rigorous
conditional infinite route, not an unconditional infinitude proof.

## 4. Consequences for the exact witness series

The repository already proves irreducibility for every odd square-family
member through `m=1057`.  Corollary 3 gives a new local proof of alternating
containment for seven of them:

| `a` | `m=2a^2-1` | `q` | `q-2` |
|---:|---:|---:|---:|
| 5 | 49 | 17 | 15 |
| 7 | 97 | 11 | 9 |
| 13 | 337 | 113 | 111 |
| 15 | 449 | 41 | 39 |
| 17 | 577 | 193 | 191 |
| 19 | 721 | 241 | 239 |
| 23 | 1057 | 353 | 351 |

The rows `a=5,13,17,19,23` belong to the cofactor-three subfamily (10).
For `m=97,337,449,577`, prime degree makes the block test especially short;
the displayed numerical criterion also verifies it directly.  The remaining
known members `m=17,161,241,881` still use the exact modular certificates.

## 5. The next local target: primes dividing `a`

Identity (1) exposes a second, presently incomplete route.  Let `p>=5`
divide `a`, put

\[
 h=v_p(L-1)=2v_p(a),\qquad L-1=up^h.
\]

Then

\[
 T(x)\equiv x(x^u-1)^{p^h}\pmod p.                  \tag{12}
\]

The exact binomial valuation

\[
 v_p\binom Lj=h-v_p(j(j-1))
 \qquad(2\le j\le p^h)                               \tag{13}
\]

shows that the Newton polygon over the residue root `1`, after division by
`(x-1)^2`, has vertices

\[
 (0,h),(p-2,h-1),(p^2-2,h-2),\ldots,(p^h-2,0).       \tag{14}
\]

Every other `u`-th root has vertices

\[
 (0,h),(p,h-1),(p^2,h-2),\ldots,(p^h,0).             \tag{15}
\]

Thus (14) contains a distinguished tame factor of degree `p-2`.  Unlike
Theorem 1, however, the competing clusters now have deeper wild layers.
Their Galois closures can have additional tame action, so the first edge
alone does **not** yet prove that a pure `(p-2)`-cycle is present.  The next
local theorem to seek is:

> prove that the tame quotient detected by the `(p-2)` edge is independent
> of the tame quotients in the later edges of (14)--(15).

Such a separation would give short structural proofs for `m=241` (`p=11`)
and `m=881` (`p=7`).  The two earlier stubborn cases `m=17` and `m=161`
have `a=3` and `9`, so this edge has length one and supplies no cycle.

## 6. Status

What is now proved is the implication

\[
 \text{irreducible}+(7)\quad\Longrightarrow\quad G=A_m,
\]

including the cofactor-three theorem (10) and seven exact witnesses.  An
unconditional infinite theorem would require either an infinite prime-value
and irreducibility result along (11), or a broader primitivity mechanism that
works without a large squarefree prime factor of `L`.  An all-irreducible
theorem additionally has to handle prime or powerful values of `L`, where
Theorem 1 supplies no proper long cycle.

The exact arithmetic identities, Newton polygons, and block inequalities in
this note are replayed by
`scripts/verify_odd_square_local_cycles.py`.
