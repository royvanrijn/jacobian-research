# Arithmetic of cancellation parameters

This is the arithmetic result attached to the cancellation construction.  It
is deliberately separate from polynomiality and from boundary geometry.
Throughout, put

\[
 n=mr,\qquad L=n+r+1,\qquad
 M_{m,r}(q)=\sum_{j=0}^{n}(-1)^j\binom Ljq^{n-j}.
\]

A root `q` selects one normalized cancellation jet.  The coefficient field
of that branch is exactly `Q(q)`: the recurrence in
[CONSTRUCTION.md](CONSTRUCTION.md) constructs every coefficient over this
field, and the constant coefficient recovers `q`.

## 1. Uniform results

### Separability and discriminant

For every `m,r>=1`,

\[
 \operatorname{disc}(M_{m,r})=
 (-1)^{n(n-1)/2}(r+1)L^{n-1}
 \binom{n+r}{r}^{n-2}.                                    \tag{1}
\]

In particular, every parameter polynomial is separable in characteristic
zero.  Formula (1) follows by passing to the affine-reciprocal polynomial

\[
 B_{n,r}(x)=\sum_{j=0}^{n}\binom{j+r}{r}x^j
\]

and using

\[
 (1-x)B'_{n,r}(x)-(r+1)B_{n,r}(x)
 =-L\binom{n+r}{r}x^n.
\]

The discriminant is a rational square precisely as follows:

- if `n=2` or `3 mod 4`, never;
- if `n=0 mod 4`, exactly when `(r+1)(n+r+1)` is a square;
- if `n=1 mod 4`, exactly when
  `(r+1) binom(n+r,r)` is a square.

For even `n`, this gives a complete parametrization.  Write `r+1=da^2` with
`d` squarefree.  Then a square discriminant occurs exactly when
`n+r+1=db^2` for some `b>a`, subject to `n=mr` and `4|n`.  Taking
`b=a+2rk` gives the infinite family

\[
 m=4dk(a+rk),\qquad k\ge1.
\]

### Proven irreducibility ranges

Irreducibility over `Q` is proved in each of the following cases:

1. `m=1`, for every `r`;
2. `n+r+1=p^k` and `v_p(n)=k-1`;
3. both `n+1` and `n+r+2` are prime and `n+r+2` is primitive modulo `n+1`;
4. `binom(n+r,r)` is prime; and
5. every pair with `mr<=30`, by exact modular factor-degree certificates.

The second item includes the useful prime case `n+r+1=p`.  The fourth follows
from the fact that every root of `B_{n,r}` lies strictly inside the unit disk.
These criteria overlap, but none is currently known to cover all pairs.

## 2. Exact finite-range Galois classification

The Galois group in its natural action is determined for every `mr<=30`.
Most cases are symmetric.  The exceptions found in that complete range are:

| Pair(s) | Group |
|---|---|
| `(2,2)` | `D_4` |
| `(2,3)`, `(6,1)` | degree-six `S_5=PGL_2(F_5)` action |
| `(4,3)` | `A_12` |
| `(2,8)`, `(16,1)` | `A_16` |
| `(17,1)` | `A_17` |
| `(1,24)`, `(12,2)` | `A_24` |

This table is an exact bounded theorem, not evidence for an all-degree
symmetric-group statement.  Formula (1) explains the alternating cases, and
also shows why alternating containment persists in infinite families.

## 3. Exact status

Proved uniformly: separability, the discriminant formula and square
criterion, the entire `m=1` irreducibility column, and the three additional
irreducibility criteria above.

Proved in a finite range: irreducibility and the complete natural Galois
group for `mr<=30`.

Open: irreducibility for all remaining `(m,r)` and an all-parameter
classification of the natural Galois group.  Minimal fields of definition of
the collision fibers are a separate arithmetic problem; see
[OPEN_PROBLEMS.md](OPEN_PROBLEMS.md).

## 4. Verification and detailed proofs

- `scripts/verify_parameter_irreducibility.py` checks the uniform criteria
  and the modular degree-sieve certificates;
- `scripts/verify_parameter_discriminant.py` checks (1), its derivation, and
  the square locus;
- `scripts/verify_parameter_galois_groups.py` and
  `scripts/verify_parameter_galois_jordan.py` certify the range `mr<=30`.

The former component proofs and tables are retained in
[the cancellation archive](../../archive/cancellation-components/).
