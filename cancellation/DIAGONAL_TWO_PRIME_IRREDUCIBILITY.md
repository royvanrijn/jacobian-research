# Two-prime Newton polygons on the divisibility diagonal

This note gives an irreducibility criterion that uses the special relation
`k=mr`.  It has two consequences not supplied by the fixed-derivative-order
theorems:

1. every fixed `m` has an explicit effective irreducible tail in `r`;
2. a general prime-window theorem gives a polynomial-scale two-parameter
   irreducibility region; and
3. the 1000 columns `1<=m<=1000` are irreducible for every `r>=1`.

Put

\[
 k=mr,\qquad N=k+r+1,
 \qquad P_{N,k}(x)=\sum_{j=0}^k\binom Njx^j.
\]

Irreducibility is unchanged by translation, so write

\[
 T_{N,k}(x)=P_{N,k}(x-1)=\sum_{j=0}^k c_jx^j.
\]

The alternating binomial identity gives

\[
 c_j=(-1)^{k-j}\binom Nj\binom{N-j-1}{k-j}
 =\frac{(-1)^{k-j}N(N-1)\cdots(N-k)}
 {j!(k-j)!(N-j)}.                                    \tag{1}
\]

This is the translated form used by
[Filaseta--Kumchev--Pasechnik](https://arxiv.org/abs/math/0409523).

## 1. A prime in the diagonal interval determines every factor degree

### Lemma

Let `p=N-u` be prime with `k<p<N`.  Then `1<=u<=r`.  If `P_(N,k)` is
reducible over `Q`, it is the product of exactly two irreducible polynomials,
of degrees

\[
 \boxed{u\quad\hbox{and}\quad k-u.}                  \tag{2}
\]

### Proof

Because `p>k`, the product of the `k+1` consecutive integers from `N-k`
through `N` contains exactly one multiple of `p`, namely `N-u=p`, and every
factor in `j!(k-j)!` is a `p`-adic unit.  Formula (1) therefore gives

\[
 v_p(c_j)=
 \begin{cases}
 0,&j=u,\\
 1,&j\ne u.
 \end{cases}                                         \tag{3}
\]

The lower Newton polygon of `T_(N,k)` consequently has the two edges

\[
 (0,1)\longrightarrow(u,0)
 \longrightarrow(k,1).                               \tag{4}
\]

Each edge has vertical displacement one, so it has no interior lattice
point.  Dumas's product theorem says that the Newton polygon of a factor is
formed from translated whole edges of (4).  Thus a proper factorization has
exactly the two degrees in (2).  Translation back from `T_(N,k)` to
`P_(N,k)` proves the claim.  QED

This lemma is also a useful conditional classification: one prime in the
interval `(k,N)` permits only one two-factor degree pattern, not an arbitrary
partition of `k`.

## 2. The two-prime criterion

### Theorem

Let `m>=2` and `r>=1`.  If the open interval

\[
 \boxed{mr<p<(m+1)r+1}                               \tag{5}
\]

contains two distinct primes, then
`P_((m+1)r+1,mr)` is irreducible over `Q`.

### Proof

Suppose the primes are `p=N-u` and `q=N-v`, where `u!=v`.  Both `u` and `v`
belong to `[1,r]`.  Since `m>=2`, this interval lies in `[1,k/2]`; hence the
unordered pairs

\[
 \{u,k-u\},\qquad \{v,k-v\}
\]

are distinct.  If the polynomial were reducible, the lemma applied at `p`
and at `q` would force both pairs to be its irreducible-factor degrees, a
contradiction.  QED

Together with the existing Eisenstein criterion when `N` is prime, this
gives the following necessary condition for a counterexample with `m>=2`:

\[
 \boxed{N\text{ is composite and }(mr,(m+1)r+1)
 \text{ contains at most one prime}.}                 \tag{6}
\]

If it contains exactly one prime `N-u`, the only possible reducibility
pattern is the two-factor pattern `(u,mr-u)`.

## 3. Every fixed-`m` column is effectively eventually irreducible

Fix `m>=2`.  The prime number theorem gives

\[
 \pi((m+1)r)-\pi(mr)\sim\frac r{\log r}
 \qquad(r\longrightarrow\infty).                    \tag{7}
\]

The interval in (5) therefore contains at least two primes for all
sufficiently large `r`.  The two-prime theorem proves

\[
 \boxed{\text{For every fixed }m,\quad
 P_{(m+1)r+1,mr}\text{ is irreducible for all sufficiently large }r.}
                                                               \tag{8}
\]

This is complementary to density-one irreducibility on every fixed-`r` row:
here the quotient `k/r=m` is fixed and the conclusion has only finitely many
exceptions.

There is a general quantitative form of this argument.

### Prime-window transfer theorem

Suppose that `H(x)>0` and that every real `x>=x_0` has a prime in
`(x,x+H(x)]`.  Put `k=mr`.  If `k>=x_0` and

\[
 H(k)+H(k+H(k))\le r,                                \tag{8a}
\]

then `P_((m+1)r+1,mr)` is irreducible.

Indeed, choose a prime `p` in `(k,k+H(k)]`, then apply the same hypothesis at
`k+H(k)` to obtain a prime `q` in
`(k+H(k),k+H(k)+H(k+H(k))]`.  The primes are distinct and (8a) puts both in
`(k,k+r+1)`.  The two-prime theorem applies.  This transfer theorem can use
any future uniform prime-gap estimate without changing the algebraic
argument.

For example, Baker--Harman--Pintz proved that, for all sufficiently large
`x`, the interval `[x-x^(21/40),x]` contains a prime; see
[The difference between consecutive primes, II](https://doi.org/10.1112/plms/83.3.532).
Apply this at `y=x+2x^(21/40)`.  Since `y-y^(21/40)>x` for all sufficiently
large `x`, every such `x` has a prime in `(x,x+2x^(21/40)]`.  Taking
`H(x)=2x^(21/40)` in (8a) gives an absolute threshold `K_0` such that

\[
 \boxed{k=mr\ge K_0,\qquad
 r\ge5k^{21/40}
 \quad\Longrightarrow\quad
 P_{(m+1)r+1,mr}\text{ is irreducible}.}             \tag{8b}
\]

Together with `k>=K_0`, the growth condition is equivalently

\[
 r\ge5^{40/19}m^{21/19}.                             \tag{8c}
\]

The threshold `K_0` is not numerically specified by the cited theorem.  This
non-explicit region is nevertheless much wider than a fixed-`m` tail: it
allows `m` to grow on the order of `k^(19/40)`.

Dusart's explicit theorem below gives a fully effective complementary
region.  Set `L=log(k)` and

\[
 \Delta(k)=\frac1{L^3}
 +\frac{1+L^{-3}}{(L+\log(1+L^{-3}))^3}.             \tag{8d}
\]

Using `H(x)=x/log^3(x)` in (8a) proves the exact effective implication

\[
 \boxed{k=mr\ge89693,\qquad m\Delta(k)\le1
 \quad\Longrightarrow\quad
 P_{(m+1)r+1,mr}\text{ is irreducible}.}             \tag{8e}
\]

Since \(\Delta(k)<2/L^3+1/L^6\), a simpler closed-form sufficient threshold is

\[
 X_m=\max\left\{89693,\,
 \exp\left((\sqrt{1+1/m}-1)^{-1/3}\right)\right\}.    \tag{8f}
\]

If `mr>X_m`, two iterations of Dusart's bound give distinct primes below
`mr(1+1/m)=(m+1)r`.  Thus (8) holds with the explicit sufficient condition

\[
 \boxed{r>X_m/m.}                                    \tag{8g}
\]

## 4. Complete columns `1<=m<=1000`

Dusart proved that for every real `a>=89693` there is a prime in

\[
 \left(a,a\left(1+\frac1{\log^3a}\right)\right].
                                                               \tag{9}
\]

See [Dusart](https://doi.org/10.1007/s11139-016-9839-4), Proposition 5.4.
Apply (9) first to `a=mr` and then to the prime it supplies.  Since the
function `1+1/log^3(a)` decreases for `a>1`, the resulting two distinct
primes are both at most

\[
 mr\rho,\qquad
 \rho=\left(1+\frac1{\log^3(89693)}\right)^2
 <1+\frac1{741}.                                    \tag{10}
\]

The strict inequality in (10) is certified with 128-bit Arb interval
arithmetic.  For every `2<=m<=741`, the last quantity is less than
`(m+1)r`; hence both
primes lie in the diagonal interval (5).  The two-prime theorem therefore
settles every pair in these columns with `mr>=89693`.

It remains to exhaust `2<=m<=741` and `mr<89693`.  First discard a pair when
`N=(m+1)r+1` is prime, by the prime-total Eisenstein criterion, or when
`(mr,N)` contains at least two primes, by Section 2.  Exactly 7689 pairs
remain: 2785 have no interval prime and 4904 have one.  Their largest degree
is 58806; the largest degree among the zero-prime cases is 34068.

Every residual pair has an exact modular factor-degree certificate.  For a
monic integral polynomial, reduction modulo any prime preserves the degree
of every rational factor; its degree must therefore be a subset sum of the
irreducible factor degrees modulo that prime, counted with multiplicity.
Intersecting these subset-sum sets over finitely many primes leaves only
`0` and `mr`.  In a one-prime case the Newton-polygon lemma makes the
certificate substantially shorter: if the interval prime is `N-u`, a
reducible polynomial has a factor of degree `u`.  At a squarefree auxiliary
prime, distinct-degree factorization only through degree `u` determines
whether the modular factors have a subset of that degree.  Here `u<=68`, so
the high-degree complementary factor never has to be factored.  The
computation through `m=300` is replayed by
`scripts/verify_parameter_irreducibility.py`; the heavier
`301<=m<=499` frontier is replayed in parallel by
`scripts/verify_parameter_irreducibility_dusart_frontier.py`, and the
`500<=m<=741` frontier by
`scripts/verify_parameter_irreducibility_sharp_dusart_frontier.py`.
The last frontier consists of 2899 residual pairs, uses auxiliary primes at
most 127, and needs at most 20 good modular steps per pair.

The already-known truncated-binomial theorem supplies `m=1`.  Consequently

\[
 \boxed{P_{(m+1)r+1,mr}(x)\text{ is irreducible for every }
 1\le m\le741,\ r\ge1.}                              \tag{11}
\]

The bound at the common threshold is uniform through `m=741`.  This is the
exact endpoint of that non-adaptive estimate: the same 128-bit Arb
calculation certifies

\[
 \rho<1+\frac1{741},
 \qquad
 \rho>1+\frac1{742}.
\]

For `742<=m<=1000`, instead begin the analytic tail at the first integer `r`
for which

\[
 \left(1+\frac1{\log^3(mr)}\right)^2<1+\frac1m.      \tag{10a}
\]

This begins at `r=122` for `m=742` and at `r=297` for `m=1000`; monotonicity
then handles every larger `r`.  Exact prime enumeration below these
column-dependent cutoffs leaves 3335 residual pairs: 1294 have no interval
prime and 2041 have one.  Their largest degree is 58823, the largest
zero-prime degree is 44298, and the largest forced degree is 54.  Exact
modular factor-degree certificates settle all of them, using auxiliary
primes at most 107 and at most 18 good steps.  They are replayed by
`scripts/verify_parameter_irreducibility_adaptive_dusart_frontier.py`.

Consequently (11) strengthens to

\[
 \boxed{P_{(m+1)r+1,mr}(x)\text{ is irreducible for every }
 1\le m\le1000,\ r\ge1.}                             \tag{11a}
\]

This is a finite-width theorem, not a proof of the full diagonal.  Formula
(8g) does, however, give a completely effective tail in every later column,
while (8b) gives a much wider asymptotic two-parameter region.
Extending the all-`r` block beyond 1000 requires only a further finite
certificate frontier, a stronger prime-interval input, or another
irreducibility criterion.

Combining (6), (8b), (8g), and (11a), any remaining counterexample must have
`m>=1001`, `r<=X_m/m`, composite `N`, and at most one prime in `(mr,N)`.
When `mr>=K_0`, it must additionally satisfy `r<5(mr)^(21/40)`.

## 5. Exact regression

Run

```bash
.venv/bin/python scripts/verify_parameter_irreducibility.py
.venv/bin/python scripts/verify_parameter_irreducibility_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_sharp_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_adaptive_dusart_frontier.py
```

The regression checks (1) and the V-shaped Newton valuations (3)--(4) on a
bounded two-prime grid.  It constructs exact modular factor-degree
certificates for the 7689 residual pairs in the common-threshold range
`2<=m<=741`, `mr<89693`, and for all 3335 residual pairs below the adaptive
tails in `742<=m<=1000`, as well as retaining the earlier complete `mr<=30`
replay.  The scripts verify the algebraic identities, finite prime
enumeration, Arb endpoint inequalities, and modular certificates.  The
adaptive replay deliberately uses the simpler upper bound
`(1+1/log^3(k))^2`; its finite frontier is therefore a conservative superset
of the frontier from the sharper effective condition (8e).  The prime number
theorem, Baker--Harman--Pintz, and Dusart are the cited analytic inputs.
