# Irreducibility results for C24 parameter polynomials

For `m,r>=1`, put

\[
 n=mr,\qquad L=n+r+1=r(m+1)+1
\]

and recall the C24 parameter polynomial

\[
 \mathcal M_{m,r}(q)=
 \sum_{j=0}^{n}(-1)^j\binom Lj q^{n-j}.                 \tag{1}
\]

This note proves irreducibility for an infinite arithmetic subfamily.  It
does not claim irreducibility for every `(m,r)`.

Separability is no longer part of this open problem: the uniform nonzero
[parameter-discriminant formula](PARAMETER_DISCRIMINANT.md) settles it for all
`m,r`.

## Complete `m=1` column from truncated-binomial theory

For positive integers `N>k`, write

\[
 P_{N,k}(x)=\sum_{j=0}^k\binom Njx^j.
\]

The C24 polynomial is its reciprocal transform:

\[
 \mathcal M_{m,r}(q)=q^{mr}P_{mr+r+1,mr}(-q^{-1}).
\]

Khanduja, Khassa, and Laishram prove that `P_(N,k)` is irreducible over `Q`
whenever

\[
 2\le 2k\le N<(k+1)^3.
\]

See [Some irreducibility results for truncated binomial
expansions](https://arxiv.org/abs/1306.0758).  For `m=1`, one has
`k=r` and `N=2r+1`, so

\[
2\le2r\le2r+1<(r+1)^3
\]

for every `r>=1`.  Reciprocal transformation preserves irreducibility;
therefore

\[
\boxed{\mathcal M_{1,r}\text{ is irreducible for every }r\ge1.}
\]

This is the first complete infinite coordinate column not restricted by a
prime-value condition.  The same numerical range also includes `(m,r)=(2,1)`,
already covered by the elementary criteria below.  For `m>=2,r>=2`, its
condition `2mr<=mr+r+1` fails, so it does not settle the remaining quadrant.

## Prime-power Eisenstein theorem

### Theorem

Suppose

\[
 L=p^k,\qquad v_p(n)=k-1                              \tag{2}
\]

for a prime `p`.  Then `M_(m,r)` is Eisenstein at `p`, hence irreducible over
`Q`.  In particular it is irreducible whenever `r(m+1)+1` is prime.

### Proof

For `0<j<p^k`, write

\[
 \binom{p^k}{j}={p^k\over j}\binom{p^k-1}{j-1}.        \tag{3}
\]

The second factor is a `p`-adic unit.  Indeed, in its product formula every
factor `(p^k-i)/i`, `1<=i<j`, is a unit, since numerator and denominator have
the same `p`-adic valuation.  Therefore

\[
 v_p\binom{p^k}{j}=k-v_p(j).                            \tag{4}
\]

Every nonleading coefficient of (1) corresponds to `1<=j<=n<L` and is
divisible by `p`.  The constant coefficient corresponds to `j=n`; condition
(2) and (4) give valuation exactly one.  The leading coefficient is one.
Eisenstein's criterion proves the theorem.  When `k=1`, condition (2) is
automatic because `0<n<p`.  QED

Thus all roots belonging to a certified pair `(m,r)` form one arithmetic
Galois orbit over `Q`.  The theorem does not decide whether two conjugate
maps become polynomially left--right equivalent over an algebraic closure;
it removes the possibility of distinct rational irreducible parameter
branches in this subfamily.

## Cyclotomic reduction theorem

### Theorem

Suppose

\[
 p=L+1\text{ is prime},\qquad
 \ell=n+1\text{ is prime},\qquad
 \operatorname{ord}_{\ell}(p)=\ell-1.                 \tag{5}
\]

Then `M_(m,r)` is irreducible over `Q`.

### Proof

Since `L=p-1` and `0<=j<=n<L`, one has

\[
 \binom Lj\equiv\binom{-1}j=(-1)^j\pmod p.
\]

Consequently

\[
 \mathcal M_{m,r}(q)\equiv
 1+q+\cdots+q^n=\Phi_\ell(q)\pmod p.                   \tag{6}
\]

Over `F_p`, the irreducible factors of the prime cyclotomic polynomial
`Phi_ell` all have degree `ord_ell(p)`.  Condition (5) therefore makes (6)
irreducible of degree `ell-1=n`.  Reduction modulo `p` preserves the degree
of the monic polynomial (1), so irreducibility over `F_p` proves
irreducibility over `Q`.  QED

This criterion is independent of the prime-power theorem.  For example,
`(m,r)=(4,1)` has `(n,L,p,ell)=(4,6,7,5)`, and `(8,5)` has
`(40,46,47,41)`; both have composite `L` and irreducible cyclotomic
reductions.

## Unit-disk leading-prime theorem

The affine-reciprocal transform of (1), up to a nonzero rational scalar, is

\[
 B_{n,r}(x)=\sum_{j=0}^n\binom{j+r}{r}x^j.              \tag{7}
\]

Indeed, if `J_(n,r)` denotes the beta-integral normalization of the parameter
polynomial, then

\[
 x^nJ_{n,r}(1-x^{-1})={n!r!\over(n+r+1)!}B_{n,r}(x).
                                                                    \tag{8}
\]

### Theorem

If `binom(n+r,r)` is prime, then `M_(m,r)` is irreducible over `Q`.  In
particular, every `r=1` member with `m+1` prime is irreducible.

### Proof

The coefficients of (7) are positive and consecutive ratios satisfy

\[
 {\binom{j+r}{r}\over\binom{j+r+1}{r}}
 ={j+1\over j+r+1}\le {n\over n+r}<1.                  \tag{9}
\]

Enestrom--Kakeya therefore puts every root of `B_(n,r)` strictly inside the
unit disk.  The polynomial is primitive and has constant coefficient one.
If it factored in `Z[x]` as two nonconstant polynomials, both factors would
have constant coefficient of absolute value one.  For either factor, the
product of the absolute values of its roots equals the absolute value of its
constant coefficient divided by its leading coefficient.  Since every root
has absolute value below one, each factor would have integer leading
coefficient of absolute value greater than one.  Their product cannot be the
prime leading coefficient `binom(n+r,r)`.  Thus (7), and hence (1), is
irreducible.  For `r=1`, the leading coefficient is `n+1=m+1`.  QED

This criterion covers cases missed by both earlier theorems, such as
`(m,r)=(6,1)`: here `m+1=7` is prime, while `L=8` fails the prime-power
valuation condition and `L+1=9` is not prime.

## Exact modular degree-sieve certificates

The three theorems above are uniform but do not cover every small pair.  The
remaining low-parameter cases admit compact exact certificates that are
stronger than asking a computer algebra system to factor over `Q`.

Let `F in Z[q]` be monic of degree `n`.  For a prime `p`, factor its reduction
as irreducibles of degrees `d_(p,1),...,d_(p,s)`, repeated according to
multiplicity, and let `S_p` be their subset-sum set.  If `F` had a monic
rational factor of degree `d`, Gauss's lemma and reduction modulo `p` would
give

\[
 d\in S_p
\]

for every `p`.  Consequently

\[
 \bigcap_{p\in\mathcal P}S_p=\{0,n\}                   \tag{10}
\]

is an exact irreducibility certificate over `Q`.

The regression constructs certificates (10) for every pair with `mr<=30`.
This range is aligned with the exact Galois classification in
[PARAMETER_GALOIS_GROUPS.md](PARAMETER_GALOIS_GROUPS.md).  It is a finite
theorem for the stated degree range, not an all-degree extrapolation.

Run the exact regression with

```bash
.venv/bin/python scripts/verify_parameter_irreducibility.py
```

It checks the coefficient valuations, cyclotomic reductions, and
affine-reciprocal transforms for certified pairs in bounded boxes, includes
every displayed C24 example, and constructs the modular degree-sieve
certificates above.  The uniform results are the three preceding proofs.
