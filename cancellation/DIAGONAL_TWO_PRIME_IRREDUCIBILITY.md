# Two-prime Newton polygons on the divisibility diagonal

This note gives an irreducibility criterion that uses the special relation
`k=mr`.  It has two consequences not supplied by the fixed-derivative-order
theorems:

1. for every fixed `m`, all sufficiently large `r` are irreducible; and
2. the six columns `1<=m<=6` are irreducible for every `r>=1`.

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

## 3. Every fixed-`m` column is eventually irreducible

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

## 4. Complete columns `1<=m<=6`

Rohrbach and Weis proved that every integer `a>=118` has a prime in

\[
 (a,14a/13);
\]

see
[Rohrbach--Weis](https://eudml.org/doc/150637), *Zum finiten Fall des
Bertrandschen Postulats*, J. Reine Angew. Math. 214/215 (1964), 432--440.
Apply their theorem first to `a=mr` and then to the prime it supplies.  This
produces two distinct primes between

\[
 mr\quad\hbox{and}\quad
 \left(\frac{14}{13}\right)^2mr=\frac{196m}{169}r.   \tag{9}
\]

For `2<=m<=6`, one has `196m<169(m+1)`.  Hence, whenever `mr>=118`, both
primes in (9) lie in the diagonal interval (5), and the two-prime theorem
applies.

The remaining finite range `2<=m<=6` and `mr<118` is proved by the exact
modular degree-sieve certificates in
`scripts/verify_parameter_irreducibility.py`.  The already-known
truncated-binomial theorem supplies `m=1`.  Consequently

\[
 \boxed{P_{(m+1)r+1,mr}(x)\text{ is irreducible for every }
 1\le m\le6,\ r\ge1.}                                \tag{10}
\]

The cutoff `6` is intrinsic only to the explicit ratio used here:
`(14/13)^2<1+1/m` holds exactly for `m<=6`.  Better explicit short-interval
results can enlarge the uniform column range, but their larger starting
points also require a correspondingly larger finite certificate computation.

## 5. Exact regression

Run

```bash
.venv/bin/python scripts/verify_parameter_irreducibility.py
```

The regression checks (1) and the V-shaped Newton valuations (3)--(4) on a
bounded two-prime grid.  It also constructs exact modular factor-degree
certificates for every pair with `mr<=30` and for every finite endpoint pair
`2<=m<=6`, `mr<118`.  The script verifies the algebraic identities and finite
certificates; the prime number theorem and the Rohrbach--Weis interval theorem
are the cited analytic inputs.
