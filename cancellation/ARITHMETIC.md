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

Filaseta--Moy also record two geometric origins of these same polynomials:
they occur in passport-one tree dessins through associated Shabat
polynomials, and the reciprocal truncated-binomial form occurs in Schubert
calculus on Grassmannians.  These connections do not currently add a uniform
arithmetic theorem, but they give two further communities standard entry
points to the family; see
[Filaseta--Moy](https://doi.org/10.4064/cm7474-3-2018), Section 1.

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

There is also a particularly simple odd-degree family not visible in (8).
For `r=1` and `k=m=1 mod 4`, the remaining square condition is

\[
 2(m+1)\text{ is a square}.
\]

Consequently

\[
 \boxed{(m,r)=(2a^2-1,1),\qquad a\text{ odd}}          \tag{8a}
\]

has square discriminant.  This begins with
`m=1,17,49,97,161,241,...` and supplies an odd-degree test family for the
alternating-group branch of the Galois question.

For fixed `r`, the even square locus is a finite union of quadratic
sequences, hence contributes `O_r(sqrt(X))` parameters `m<=X`.  If `r` is
odd and at least three, the odd square condition is the integral-point
equation

\[
 Y^2=(r+1)\binom{r(m+1)}r,                             \tag{8b}
\]

whose squarefree right side has degree `r`; the resulting hyperelliptic
curve has genus `(r-1)/2`.  Siegel's theorem makes this odd branch finite.
For even `r` there is no odd-degree branch, while `r=1` is exactly (8a).
Therefore

\[
 \#\{m\le X:\operatorname{disc}(M_{m,r})\text{ is a square}\}
 =O_r(\sqrt X).                                        \tag{8c}
\]

Square discriminants are thus infinite but density zero on every fixed-`r`
row.  A detailed derivation is retained in
[the discriminant note](../archive/cancellation-components/PARAMETER_DISCRIMINANT.md).

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

The derivative form (3a) nevertheless gives a strong asymptotic theorem on
every fixed-`r` row.  Borisov--Filaseta--Lam--Trifonov prove that, for each
fixed derivative order `r`, all but

\[
 O\!\left(\frac{T\log\log T}{\log T}\right)
\]

geometric degrees `h<=T` give an irreducible `G_h^(r)`.  Restricting their
exceptional set to `h=r(m+1)` yields

\[
 \boxed{\#\{m\le X:M_{m,r}\text{ reducible}\}
 =O_r\!\left(\frac{X\log\log X}{\log X}\right).}       \tag{10a}
\]

In particular, for every fixed `r`, irreducibility holds for a density-one
set of `m`.  For `r=1`, their sharper first-derivative theorem replaces the
right side by `O_(epsilon)(X^(1/3+epsilon))`.  See
[Borisov--Filaseta--Lam--Trifonov](https://doi.org/10.4064/aa-90-2-121-153).

The same paper records additional proved first-derivative families.  Since
the geometric degree is `h=m+1` when `r=1`, irreducibility is known if

\[
 h=p-1,\qquad h=p^a,\qquad h+1\text{ is squarefree},
 \qquad\text{or}\qquad h=2p-1,                         \tag{10b}
\]

where `p` is prime and `a>=1`.  The first case is the existing prime case
`N=p`; the other three enlarge the explicitly proved `r=1` subfamilies.

The Jacobi form (3b) supplies a complementary fixed-degree theorem.
Cullinan--Hajir--Sell prove that for every fixed `k>=6`, `J_k(x,y_0)` is
irreducible for all but finitely many rational `y_0`.  For odd `k`, its
Galois group is `S_k` outside a finite set; for even `k`, the `A_k`
specializations lie in a thin set.  They also prove that the generic
group over `Q(y)` is `S_k`.  See
[Cullinan--Hajir--Sell](https://doi.org/10.5802/jtnb.659).  This does not by
itself settle the moving specialization `(k,y)=(mr,r)`, but it identifies the
relevant specialization geometry and explains why large Galois groups are
the generic expectation.

Combining (8c) and (10a) shows that, for every fixed `r`, a density-one set
of `m` gives an irreducible polynomial with nonsquare discriminant.  Hence a
density-one proof that the Galois group contains `A_k` would automatically
upgrade to the sharper conclusion `Gal(P_(N,k))=S_k` for density-one `m`.

The proof of the geometric-derivative theorem contains one further uniform
Galois consequence.  Extracting its `p`-adic clusters shows that, for fixed
`r`, outside

\[
 O_r\!\left(\frac{X\log\log X}{\log X}\right)
\]

parameters `m<=X`, the transitive Galois group has order divisible by a prime

\[
 (\log X)^{10}<p\le mr-3.                              \tag{10c}
\]

Equivalently, it contains an order-`p` element of cycle type `p^t 1^f` with
`t>=1` and `f>=p-1`.  This is genuine growing wild ramification, but it is not
yet a Jordan-cycle theorem: several root-of-unity residue clusters may be
ramified simultaneously, so the argument does not prove `t=1`.  The exact
extraction and this support warning are proved in
[FIXED_R_NEWTON_RAMIFICATION.md](FIXED_R_NEWTON_RAMIFICATION.md).

The divisibility diagonal also supports a complementary Newton-polygon
theorem in the other direction.  If `m>=2` and the interval

\[
 mr<p<(m+1)r+1
\]

contains two primes, the two corresponding translated Newton polygons force
two different and incompatible factor-degree pairs.  Hence the polynomial is
irreducible.  The prime number theorem shows that, for every fixed `m`, this
criterion holds for all sufficiently large `r`.  An explicit prime-interval
theorem of Rohrbach--Weis, combined with finite exact certificates below
degree `118`, proves the complete columns

\[
 \boxed{1\le m\le6,\qquad r\ge1.}                    \tag{10d}
\]

More precisely, if a counterexample with `m>=2` exists, then `N` is composite
and `(mr,(m+1)r+1)` contains at most one prime.  If it contains exactly one
prime `N-u`, the polynomial must be the product of exactly two irreducibles
of degrees `u` and `mr-u`.  The proof and the short-interval calculation are
given in
[DIAGONAL_TWO_PRIME_IRREDUCIBILITY.md](DIAGONAL_TWO_PRIME_IRREDUCIBILITY.md).

In addition, the repository proves irreducibility on the cancellation
diagonal in each of the following cases:

1. `m=1`, for every `r`, by the first cited truncated-binomial theorem;
2. `N=p^a` and `v_p(k)=a-1`, by an explicit Eisenstein argument;
3. both `k+1` and `N+1` are prime and `N+1` is primitive modulo `k+1`, by
   cyclotomic reduction;
4. `binom(N-1,k)` is prime, using the strict unit-disk location of the roots
   of `B_(k,r)`;
5. the interval `(mr,(m+1)r+1)` contains at least two primes, by incompatible
   translated Newton polygons; and
6. every pair with `mr<=30`, together with every finite endpoint
   `2<=m<=6`, `mr<118`, by exact modular factor-degree certificates.

The second item includes the useful prime case `N=p`.  The bounded
certificates in item 6 are independently replayable even where their
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

For a fixed omitted-term parameter `r`, (10c) now supplies large prime-order
inertia on a density-one set.  What remains is sharper and more local than
the earlier formulation suggested: isolate one `p`-cycle from the possibly
simultaneous `p`-clusters, then rule out nontrivial block systems.  Only after
both steps does Jordan's theorem yield `A_(mr)`.  In particular, a Newton
edge of length `p` must not be described as a pure `p`-cycle without a support
argument.

The natural Galois group is determined in this repository for every
`mr<=30`, and also for two further members of the odd square family.  Most
cases in the complete degree range are symmetric.  The known non-symmetric
groups are:

| Pair(s) | Group |
|---|---|
| `(2,2)` | `D_4` |
| `(2,3)`, `(6,1)` | degree-six `S_5=PGL_2(F_5)` action |
| `(4,3)` | `A_12` |
| `(2,8)`, `(16,1)` | `A_16` |
| `(17,1)` | `A_17` |
| `(1,24)`, `(12,2)` | `A_24` |
| `(49,1)` | `A_49` |
| `(97,1)` | `A_97` |

The two degree-six exceptions reproduce the published Filaseta--Moy
phenomenon.  Formula (5) explains every alternating entry and the infinite
square families.  The degree-49 and degree-97 entries are the next two
members of the odd family (8a), proved by exact modular degree sieves and
Jordan-cycle certificates.  The data therefore suggest the sharpened
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
the six complete columns `1<=m<=6`, and the additional diagonal
irreducibility criteria above.  The geometric-derivative literature further
establishes density-one irreducibility on every fixed-`r` row, with the
quantitative bounds in (10a), and the extra `r=1` families (10b).  Exact
finite certificates establish irreducibility for `mr<=30` and for the finite
endpoints `2<=m<=6`, `mr<118`.  The complete natural Galois group is known
for `mr<=30`, together with `A_49` and `A_97` on the odd square family.

Open: prove (10) in the remaining cases and answer the large-group question
(12).  Minimal fields of definition of the collision fibers are a separate
arithmetic problem; see the [research roadmap](RESEARCH_ROADMAP.md).

The exact regressions are:

- `scripts/verify_parameter_irreducibility.py`, for the uniform diagonal
  criteria, two-prime Newton polygons, and modular degree-sieve certificates;
- `scripts/verify_parameter_discriminant.py`, for (3a), (5)--(8a) after
  specialization;
- `scripts/verify_fixed_r_newton_ramification.py`, for the reciprocal
  numerator, local congruence, and representative cyclotomic Newton edges;
- `scripts/verify_parameter_galois_groups.py` and
  `scripts/verify_parameter_galois_jordan.py`, for the range `mr<=30` and
  the additional odd square-family groups `A_49` and `A_97`.

Detailed component proofs and tables are retained in
[the cancellation archive](../archive/cancellation-components/).
