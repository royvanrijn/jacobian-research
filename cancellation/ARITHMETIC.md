# Truncated-binomial arithmetic of cancellation parameters

This is the arithmetic result attached to the cancellation construction.  It
is deliberately separate from polynomiality and from boundary geometry.  To
place it in standard number-theoretic notation, write

\[
 P_{N,k}(x)=\sum_{j=0}^{k}\binom Njx^j,
 \qquad
 N=(m+1)r+1,
 \qquad
 k=mr=N-r-1.
\]

Then the cancellation parameter polynomial is exactly the signed reciprocal
transform

\[
 \boxed{M_{m,r}(q)=q^kP_{N,k}(-q^{-1}).}                 \tag{1}
\]

Thus reciprocity identifies irreducibility, root fields, splitting fields,
and the natural permutation Galois groups of `M_(m,r)` and `P_(N,k)`.  The
parameter family is the divisibility diagonal

\[
 1\le k\le N-2,
 \qquad
 N-k-1\mid k.                                           \tag{2}
\]

Equivalently, it retains all terms through degree `N-r-1` and omits the last
`r+1` terms of the full binomial expansion.

The affine-reciprocal polynomial already used in this repository is also a
standard form:

\[
 B_{k,r}(x)=\sum_{j=0}^{k}\binom{r+j}{r}x^j,
 \qquad
 x^kM_{m,r}(1-x^{-1})=(-1)^kB_{k,r}(x).                 \tag{3}
\]

It is precisely Filaseta--Moy's polynomial `p_(d,t)` with `(d,t)=(k,r)`.
There are two further classical descriptions.  If

\[
 G_h(x)=1+x+\cdots+x^h,
 \qquad h=N-1=(m+1)r,
\]

then termwise differentiation gives the exact identity

\[
 \boxed{B_{k,r}(x)=\frac1{r!}G_{N-1}^{(r)}(x).}        \tag{3a}
\]

Thus the cancellation problem is also the divisor-diagonal subfamily
`h=(m+1)r` of the classical conjecture that every nontrivial derivative of a
geometric polynomial is irreducible.  Moreover, in the one-parameter Jacobi
family

\[
 J_d(x,y)=\sum_{j=0}^d\binom{y+j}{j}x^j
\]

one has simply

\[
 B_{k,r}(x)=J_k(x,r).                                  \tag{3b}
\]

The four notational systems are therefore related by

| Cancellation notation | Truncated binomial | Geometric derivative | Jacobi / Filaseta--Moy |
|---|---|---|---|
| `L=(m+1)r+1` | outer exponent `N` | geometric degree `h=N-1` | `d+t+1=N` |
| `n=mr` | truncation degree `k` | derivative degree `h-r=k` | Jacobi degree `d=k` |
| `r` | omitted count minus one | derivative order `r` | `y=t=r` |

The hypergeometric form used elsewhere in the cancellation calculation is
equivalently the Jacobi-polynomial identity

\[
 {}_2F_1(-k,1;r+2;q)
 =\frac{k!}{(r+2)_k}P_k^{(r+1,-N)}(1-2q).             \tag{3c}
\]

The parameter `-N` moves with the degree and lies outside the usual
orthogonality range, so (3c) is structural rather than an immediate
irreducibility theorem.

A root `q` selects one normalized cancellation jet.  Its coefficient field
is exactly `Q(q)`: the recurrence in
[CONSTRUCTION.md](CONSTRUCTION.md) constructs every coefficient over this
field, and the constant coefficient recovers `q`.

## 1. Separability and discriminant

The classical identity

\[
 P_{N,k}(x)-\frac{x+1}{N}P'_{N,k}(x)
 =\binom{N-1}{k}x^k                                    \tag{4}
\]

shows immediately that `P_(N,k)` has no common root with its derivative.
This standard separability argument appears in
[Filaseta--Kumchev--Pasechnik](https://arxiv.org/abs/math/0409523).

For `k>=2`, the general truncated-binomial discriminant formula is

\[
 \operatorname{disc}(P_{N,k})=
 (-1)^{k(k-1)/2}(N-k)N^{k-1}
 \binom{N-1}{k}^{k-2}.                                  \tag{5}
\]

It is Lemma 3 of
[Filaseta--Moy](https://doi.org/10.4064/cm7474-3-2018), stated there in the
equivalent `p_(d,t)` coordinates.  The degree-one case is trivial.  Reciprocal
transformation preserves this discriminant because `P_(N,k)` has constant
coefficient one, so substituting `N=k+r+1` in (5) gives

\[
 \operatorname{disc}(M_{m,r})=
 (-1)^{k(k-1)/2}(r+1)N^{k-1}
 \binom{k+r}{r}^{k-2}.                                  \tag{6}
\]

The repository retains an independent exact derivation through (3) and

\[
 (1-x)B'_{k,r}(x)-(r+1)B_{k,r}(x)
 =-N\binom{k+r}{r}x^k.                                  \tag{7}
\]

Formula (5) makes the square criterion transparent:

- if `k=2` or `3 mod 4`, the discriminant is never a rational square;
- if `k=0 mod 4`, it is a square exactly when `N(N-k)` is a square;
- if `k=1 mod 4`, it is a square exactly when
  `(N-k) binom(N-1,k)` is a square.

On the cancellation diagonal this is the criterion formerly written in
`(n,r,L)` notation.  For even `k`, write `r+1=da^2` with `d` squarefree.
Then a square discriminant occurs exactly when `N=db^2` for some `b>a`,
subject to `k=mr` and `4|k`.  Taking `b=a+2rs` gives the infinite family

\[
 m=4ds(a+rs),\qquad s\ge1.                              \tag{8}
\]

Thus alternating containment is an infinite phenomenon on this diagonal,
not a sporadic feature of the finite table below.

## 2. Irreducibility: literature and proved diagonal ranges

The general conjecture for truncated binomial polynomials is

\[
 P_{N,k}(x)\text{ is irreducible over }\mathbb Q
 \quad(1\le k\le N-2).                                  \tag{9}
\]

It is stated in this form, together with recent progress, by
[Laishram--Yadav](https://doi.org/10.1142/S1793042124500817).  The
cancellation problem is the structured subfamily (2):

\[
 \boxed{P_{(m+1)r+1,mr}(x)\text{ is irreducible for all }m,r\ge1.} \tag{10}
\]

Known general theorems intersect this diagonal but do not cover it:

- Khanduja--Khassa--Laishram prove irreducibility when
  `2<=2k<=N<(k+1)^3`.  Here `2k<=N` becomes `(m-1)r<=1`, giving the full
  `m=1` column and the isolated pair `(2,1)`; see
  [their theorem](https://arxiv.org/abs/1306.0758).
- Dubickas--Siurys prove irreducibility for every `k<=6` and `N>=k+2`; see
  [their small-degree theorem](https://doi.org/10.1007/s12044-016-0325-0).
- The fixed-`k`, sufficiently-large-`N` theorems do not settle (10), because
  `k=mr` grows with `N` along every infinite fixed-`m` or fixed-`r` row.

In addition, the repository proves irreducibility on the cancellation
diagonal in each of the following cases:

1. `m=1`, for every `r`, by the first cited truncated-binomial theorem;
2. `N=p^a` and `v_p(k)=a-1`, by an explicit Eisenstein argument;
3. both `k+1` and `N+1` are prime and `N+1` is primitive modulo `k+1`, by
   cyclotomic reduction;
4. `binom(N-1,k)` is prime, using the strict unit-disk location of the roots
   of `B_(k,r)`; and
5. every pair with `mr<=30`, by exact modular factor-degree certificates.

The second item includes the useful prime case `N=p`.  The bounded
certificates in item 5 are independently replayable even where their
conclusions overlap published small-degree or earlier computational ranges.
None of the cited theorems or these criteria covers every pair in (10).

## 3. Galois groups and transferable methods

Filaseta--Moy prove that, for fixed truncation degree `k!=6` and sufficiently
large `N`, the Galois group of `P_(N,k)` is `S_k`; they isolate a possible
`PGL_2(F_5)` phenomenon in degree six.  Their explicit exceptional
polynomials `p_(6,1)` and `p_(6,3)` are exactly the cancellation pairs
`(m,r)=(6,1)` and `(2,3)`.

Their Newton-polygon/Jordan-cycle lemmas apply verbatim to (1)--(3), but the
global existence argument does not transfer automatically.  It holds `k`
fixed while the complementary parameter grows; on the cancellation diagonal
`k=mr` and `r` generally grow together.  For example, their large-cycle step
asks for a prime cycle length `ell` with `k/2<ell<k-2` and suitable large
prime factors of numbers that become

\[
 N-\ell,\qquad r+1+\ell.                                \tag{11}
\]

Uniformly producing such factors is a moving prime-factor problem rather
than the fixed-degree Thue-equation problem treated in their paper.

The natural Galois group is determined in this repository for every
`mr<=30`.  Most cases are symmetric.  The exceptions are:

| Pair(s) | Group |
|---|---|
| `(2,2)` | `D_4` |
| `(2,3)`, `(6,1)` | degree-six `S_5=PGL_2(F_5)` action |
| `(4,3)` | `A_12` |
| `(2,8)`, `(16,1)` | `A_16` |
| `(17,1)` | `A_17` |
| `(1,24)`, `(12,2)` | `A_24` |

The two degree-six exceptions reproduce the published Filaseta--Moy
phenomenon.  Formula (5) explains every alternating entry and the infinite
square families.  The finite data therefore suggest the sharpened
all-parameter question: apart from the displayed low-degree `D_4` and
`PGL_2(F_5)` exceptions, does irreducibility imply

\[
 \operatorname{Gal}(P_{N,k}/\mathbb Q)=
 \begin{cases}
 A_k,&\operatorname{disc}(P_{N,k})\text{ is a square},\\
 S_k,&\operatorname{disc}(P_{N,k})\text{ is not a square}
 \end{cases}                                             \tag{12}
\]

on the divisibility diagonal (2)?

## 4. Known boundary and verification

Uniform results establish separability, the specialized square criterion,
the full `m=1` irreducibility column, and the three additional diagonal
irreducibility criteria above.  Exact finite certificates establish
irreducibility and the complete natural Galois group for `mr<=30`.

Open: prove (10) in the remaining cases and answer the large-group question
(12).  Minimal fields of definition of the collision fibers are a separate
arithmetic problem; see [OPEN_PROBLEMS.md](OPEN_PROBLEMS.md).

The exact regressions are:

- `scripts/verify_parameter_irreducibility.py`, for the uniform diagonal
  criteria and modular degree-sieve certificates;
- `scripts/verify_parameter_discriminant.py`, for (5)--(8) after
  specialization;
- `scripts/verify_parameter_galois_groups.py` and
  `scripts/verify_parameter_galois_jordan.py`, for the range `mr<=30`.

Detailed component proofs and tables are retained in
[the cancellation archive](../archive/cancellation-components/).
