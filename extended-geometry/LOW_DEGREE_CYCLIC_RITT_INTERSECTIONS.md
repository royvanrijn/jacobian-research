# Low-degree cyclic symmetry--Ritt intersections

Work over an algebraically closed field of characteristic zero.  This note
computes the intersection proposed in
[INCIDENCE_COVER_AUTOMORPHISMS.md](INCIDENCE_COVER_AUTOMORPHISMS.md) for

\[
 N=6,8,9,10,12.
\]

It treats every marked cyclic stratum, not only the exact-double locus.  The
bounded exact checker is
[`verify_low_degree_cyclic_ritt_intersections.py`](../scripts/verify_low_degree_cyclic_ritt_intersections.py).

## 1. Conventions and elimination

Let `S_(e,q)` denote the exact cyclic stratum

\[
 H(W)=\sum_{k=0}^m h_kW^{e+kq},\qquad m={N-e\over q},       \tag{1.1}
\]

where `q>1`, `h_0h_m!=0`, and the affine Hessian stabilizer is exactly
`mu_q`.  The endpoint equations are

\[
 \sum h_k=0,\qquad q\sum k h_k=-1,                         \tag{1.2}
\]

and the admissible open also removes `H''(1)=-2`.  If
`I={k:h_k!=0}`, exactness is simply

\[
 \gcd(I)=1.                                                \tag{1.3}
\]

Write `D_(a,b)` for the locus where `H-s_0W=A o B`, with outer degree `a`,
inner degree `b`, and `ab=N`.  The computation normalizes `B` to be monic
with `B(0)=0`.  Its coefficients are recovered successively from degrees
`N-1` down to `N-b+1`; the remaining coefficients are the decomposition
equations.  For every exact support (1.3), the checker saturates by

\[
 h_0h_m\bigl(H''(1)+2\bigr)\prod_{k\in I}h_k.              \tag{1.4}
\]

Consequently a dash in the tables is an exact empty-open certificate, not a
failure to find a small rational point.

## 2. Complete census

At the level of the **union** of cyclic strata, the answer is immediate and
uniform: every ordered vertical Ritt locus in the five degrees meets.  If
`N=ab`, then

\[
 \boxed{H_{N,b}(W)={W^b-W^N\over N-b}}
 \in S_{b,N-b}\cap D_{a,b},                         \tag{2.1}
\]

because

\[
 H_{N,b}(W)={Y-Y^a\over N-b}\bigg|_{Y=W^b}.         \tag{2.2}
\]

This seed has `H''(1)=-(N+b-1)!=-2`, and its exact cyclic stabilizer is
`mu_(N-b)`.  Thus no factor pair in `N=6,8,9,10,12` is disjoint from the
cyclic Hessian-symmetry locus.  The nontrivial refinement is to determine
which **individual exact** strata it meets; that is the census below.

In every table, columns are ordered by outer degree first.  A star marks a
nonempty cell not explained by the direct monomial condition
`b|gcd(e,q)`.

### Degree 6

| `(e,q)` | `D_(3,2)` | `D_(2,3)` |
|---|---:|---:|
| `(2,2)` | yes | -- |
| `(2,4)` | yes | -- |
| `(3,3)` | -- | yes |
| `(4,2)` | yes | -- |

### Degree 8

| `(e,q)` | `D_(4,2)` | `D_(2,4)` |
|---|---:|---:|
| `(2,2)` | yes | yes* |
| `(2,3)` | -- | -- |
| `(2,6)` | yes | -- |
| `(3,5)` | -- | -- |
| `(4,2)` | yes | -- |
| `(4,4)` | yes | yes |
| `(5,3)` | -- | -- |
| `(6,2)` | yes | -- |

### Degree 9

| `(e,q)` | `D_(3,3)` |
|---|---:|
| `(2,7)` | -- |
| `(3,2)` | yes* |
| `(3,3)` | yes |
| `(3,6)` | yes |
| `(4,5)` | -- |
| `(5,2)` | yes* |
| `(5,4)` | -- |
| `(6,3)` | yes |
| `(7,2)` | -- |

### Degree 10

| `(e,q)` | `D_(5,2)` | `D_(2,5)` |
|---|---:|---:|
| `(2,2)` | yes | -- |
| `(2,4)` | yes | -- |
| `(2,8)` | yes | -- |
| `(3,7)` | -- | -- |
| `(4,2)` | yes | -- |
| `(4,3)` | -- | -- |
| `(4,6)` | yes | -- |
| `(5,5)` | -- | yes |
| `(6,2)` | yes | -- |
| `(6,4)` | yes | -- |
| `(7,3)` | -- | -- |
| `(8,2)` | yes | -- |

### Degree 12

| `(e,q)` | `D_(6,2)` | `D_(4,3)` | `D_(3,4)` | `D_(2,6)` |
|---|---:|---:|---:|---:|
| `(2,2)` | yes | yes* | yes* | yes* |
| `(2,5)` | -- | -- | -- | -- |
| `(2,10)` | yes | -- | -- | -- |
| `(3,3)` | -- | yes | -- | yes* |
| `(3,9)` | -- | yes | -- | -- |
| `(4,2)` | yes | -- | yes* | yes* |
| `(4,4)` | yes | -- | yes | -- |
| `(4,8)` | yes | -- | yes | -- |
| `(5,7)` | -- | -- | -- | -- |
| `(6,2)` | yes | -- | -- | -- |
| `(6,3)` | -- | yes | -- | -- |
| `(6,6)` | yes | yes | -- | yes |
| `(7,5)` | -- | -- | -- | -- |
| `(8,2)` | yes | -- | -- | -- |
| `(8,4)` | yes | -- | yes | -- |
| `(9,3)` | -- | yes | -- | -- |
| `(10,2)` | yes | -- | -- | -- |

Thus the bounded classification is

\[
 \boxed{D_{a,b}\cap S_{e,q}\ne\varnothing
 \iff b\mid\gcd(e,q),}
\]

apart from exactly nine starred cells:

\[
\begin{gathered}
 (8;2,2;4),\\
 (9;3,2;3),\ (9;5,2;3),\\
 (12;2,2;3),\ (12;2,2;4),\ (12;2,2;6),\\
 (12;3,3;6),\\
 (12;4,2;4),\ (12;4,2;6).
\end{gathered}                                             \tag{2.3}
\]

Here a tuple is `(N;e,q;b)`.

## 3. Exact-double subtable

The orbifold strata used on the exact-double marked-seed open have `e=2`.
Extracting only those rows gives the particularly clean answer:

| `N` | exact `q` | vertical Ritt loci met |
|---:|---:|---|
| 6 | 2 | `D_(3,2)` |
| 6 | 4 | `D_(3,2)` |
| 8 | 2 | `D_(4,2)`, `D_(2,4)` |
| 8 | 3 | none |
| 8 | 6 | `D_(4,2)` |
| 9 | 7 | none |
| 10 | 2, 4, 8 | `D_(5,2)` |
| 12 | 2 | all four loci |
| 12 | 5 | none |
| 12 | 10 | `D_(6,2)` |

The odd containing-symmetry rows `(N,q)=(8,3),(12,5)` do meet the quadratic-
inner locus in their closures, but only after the symmetry jumps to `q=6`
and `q=10`; their exact strata do not meet it.

## 4. Witnesses and collision types

### Direct monomial intersections

If `b|gcd(e,q)`, every member of `S_(e,q)` is a polynomial in `W^b`.
Taking `s_0=0` gives the direct inner factor `B=W^b`.  These are all of the
unstarred positive cells.

### Chebyshev intersections

For `N` divisible by four,

\[
 H_N(W)={1-T_N(W)\over N^2}                         \tag{4.1}
\]

is an admissible exact `(e,q)=(2,2)` seed.  Since
`T_N=T_a o T_b`, it belongs to every `D_(a,b)`.  Formula (4.1) supplies the
starred degree-eight cell and the three starred cells in the first
degree-twelve row.

### The two nonzero-vertical degree-nine witnesses

For `S_(3,2)` take

\[
 H={-W^9-3W^7-3W^5+7W^3\over24},\qquad s_0=-{1\over3}.
\]

Then

\[
 H-s_0W=-{(W^3+W)^3-8(W^3+W)\over24}.              \tag{4.2}
\]

For `S_(5,2)`, let `rho^2=-3` and take

\[
 H={(1-\rho)W^5+2\rho W^7-(1+\rho)W^9\over4},
 \qquad s_0={1\over18}.                             \tag{4.3}
\]

With

\[
 B=W^3-{3+\rho\over6}W,qquad
 A(Y)=-{1+\rho\over4}Y^3+{3-\rho\over36}Y,
\]

equation (4.3) satisfies `H-s_0W=A(B)`.  These two cells show that a search
restricted to `s_0=0` would give the wrong degree-nine answer.

### The remaining nested degree-twelve witnesses

The `S_(3,3) intersect D_(2,6)` cell has

\[
 H=-{W^3\over2}+{41W^6\over12}-5W^9+{25W^{12}\over12}
   ={25\over12}B^2+{5\over12}B,
 \quad B={W^3(5W^3-6)\over5}.                      \tag{4.4}
\]

The two `S_(4,2)` cells are

\[
 -{(W^4+W^2)^2(W^4+W^2-2)\over24}\in D_{3,4},     \tag{4.5}
\]

and

\[
 -{(W^6+W^4)(W^6+W^4-2)\over20}\in D_{2,6}.       \tag{4.6}
\]

All witnesses (4.1), (4.4)--(4.6) use `s_0=0`.

## 5. Uniform indication

The census separates three mechanisms cleanly:

1. the generic monomial mechanism `b|gcd(e,q)`;
2. Chebyshev collisions on the exact-double `q=2` rows when `4|N`;
3. lower-dimensional nested-composition incidences, plus the two genuinely
   vertical degree-nine points/families.

There are no additional cells through degree twelve.  In particular, odd
cyclic order by itself does not force vertical imprimitivity, and the
maximal-symmetry binomial seeds meet only those loci visible from their
actual common exponent divisor.

## Reproduction

The checker requires SymPy and Singular:

```bash
.venv/bin/python scripts/verify_low_degree_cyclic_ritt_intersections.py
```

It prints the five tables and ends with three `PASS` lines.
