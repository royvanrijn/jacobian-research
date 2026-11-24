# Split-prime sumsets as the prime discriminant varies

For a rational prime \(p\equiv3\bmod4\), let \(C_p\) consist of 1 and
the positive integers supported entirely on primes split in
\(\mathbb Q(\sqrt{-p})\). For \(p\ge7\), the admissible set \(L_p\)
consists of \(n\ge2\) with \(p\nmid n\), requiring \(n\) even if
\(p\equiv3\bmod8\). Put \(L_3=\{n\ge2:n\equiv2\bmod6\}\).
Write \(E_p=L_p\setminus(C_p+C_p)\), and \(M(p)=\max(E_p\cup\{0\})\)
when \(E_p\) is finite.

## The law of gaps in a fixed window

**Theorem.** Conditional on either \(p\equiv3\bmod8\) or
\(p\equiv7\bmod8\), the splitting indicators of any fixed finite
collection of odd primes are independent fair Bernoulli variables
in the limit over prime discriminants. The law of \(E_p\cap[2,B]\),
for each fixed integer \(B\ge4\), is the pushforward of this law for
the odd primes at most \(B-1\), with the status of 2 fixed by the
chosen progression. Among primes \(p\equiv3\bmod8\), for even
\(B\ge4\), the limiting proportion with every even \(n\) from 4
through \(B\) in \(E_p\) is exactly \(2^{1-\pi(B-1)}\).

Indeed, for odd \(q\ne p\), quadratic reciprocity gives
\[
 \left(\frac{-p}{q}\right)=\left(\frac q p\right)
                         =\left(\frac p q\right)\left(\frac{-1}q\right)
 \quad (p\equiv3\bmod4).
\]
Fixing a splitting sign therefore selects exactly \((q-1)/2\)
reduced residues for \(p\) modulo \(q\), whatever the fixed sign
\(({-1}/q)\). The Chinese remainder theorem multiplies the counts
over \(q\), even after fixing \(p\bmod8\). The prime number theorem
in arithmetic progressions gives the stated natural proportions for
each of these finitely many reduced classes. In particular every
prescribed pattern occurs infinitely often. On the larger progression
\(p\equiv3\bmod4\), the two possible statuses of 2 are equally likely
and independent of the odd-prime signs in this finite sense.
[Soprounov, Theorem 2.1, fixed-modulus prime number theorem](../foundations/extended-geometry/sources/soprounov_pnt_arithmetic_progressions.pdf)

For \(p>B\), the admissibility test in the window depends only on the
fixed status of 2. A summand of a positive representation of
\(n\le B\) is at most \(B-1\). Thus only the splitting bits for primes
up to \(B-1\) matter. The limiting probability of any specified gap
set in the window is the number of odd-prime bit vectors producing it
divided by \(2^{\pi(B-1)-1}\).

When 2 is inert, if all odd primes below \(B\) are inert, the only
core value at most \(B-1\) is 1. Conversely a split prime \(q<B\)
gives \(q+1=1+q\), an even represented target between 4 and \(B\).
This proves the run formula and gives arbitrarily long initial runs
of eligible gaps as \(p\) varies. If 2 splits, \(4=2+2\).

The limiting gap indicators are positively associated. A gap
indicator decreases when more splitting bits are turned on, so two
increasing functions of the gap vector are both decreasing functions
of the independent bits. Their covariance is nonnegative by induction
on the number of bits. Condition on the last bit to apply the
induction hypothesis to the conditional covariances. The two
conditional means are decreasing functions of that same bit, so
their covariance is also nonnegative. For example, on
\(p\equiv3\bmod8\) the limiting gap probabilities at 4 and 6 are
\(1/2\) and \(1/4\), and their joint gap probability is \(1/4\).
The bit for 3 controls 4, and 6 is a gap precisely when 3 and 5
are both inert. These finite-window laws concern fixed \(B\) as
\(p\) grows.

## An effective polynomial bound

**Theorem.** For every prime \(p\equiv3\bmod4\), the set \(E_p\)
is finite and
\[
                    M(p)< 2^{3696}p^4 .
\]
More generally, for each integer \(h>4\), put
\[
 T_h=\left(
   \prod_{\substack{q<2^h\\q\text{ prime}}}
   \max_{0\le a\le2h}\frac{(a+1)^h}{q^a}
   \right)^{1/h}.
\]
Every \(n\in L_p\) with
\[
                n\ \ge\ (256T_h)^{4h/(h-4)}p^{3h/(h-4)}
                                                              \tag{1}
\]
is represented. In particular \(M(p)\ll_\varepsilon p^{3+\varepsilon}\)
effectively for every fixed \(\varepsilon>0\). No class number
hypothesis is imposed.

The admissibility restrictions are necessary. All core values are
nonzero squares modulo \(p\), and \(-1\) is a nonsquare. In the inert
dyadic case the core values are odd. Conversely for \(p\ge7\) and
\(p\nmid n\), there are \(p+1-2(1+(n/p))>0\) ordered solutions of
\(u^2+v^2=n\bmod p\) with both coordinates nonzero. This follows by
counting the \(p+1\) elements in each nonzero norm fiber of
\(\mathbb F_{p^2}/\mathbb F_p\) and removing the two axes.
At any other odd prime one can choose both summands to be units,
since this excludes at most two choices. These choices lift and may
be realized by the Chinese remainder theorem and split primes in
reduced progressions. When 2 splits, admit an even summand for an odd
target by taking it twice an odd split prime. The factor 2 preserves
the nonzero square classes modulo \(p\), so the same residue argument
applies. At \(p=3\) the nonzero square set modulo 3 is \(\{1\}\),
giving \(L_3\). The exact earlier identities give \(E_3=\varnothing\)
and \(E_7=\{11,29\}\), which meet (1).
[Eisenstein core identity](../../results/eisenstein-restricted-prime-sums/sections/primitive-sumsets.tex) ·
[Discriminant −7](PRIMITIVE_NORM_MINUS7_EXACT.md)

### A shell estimate with the discriminant free

Let \(p\ge11\) and \(p\equiv3\bmod4\). Set
\(R=\mathbb Z[(1+\sqrt{-p})/2]\) and
\[
 Q_p(x,y)=x^2+xy+(p+1)y^2/4,\qquad O=R\oplus Rj,\qquad
 j^2=-1,\quad jz=\overline z j .
\]
This is a maximal order in the definite quaternion algebra of
discriminant \(p\). Its norm is \(Q_p\oplus Q_p\). The norm-one
units are \(1,-1,j,-j\), since the only field norm-one elements
are \(1,-1\). The inverse-projective-unit mass is \((p-1)/12\)
and the probability of the principal class is therefore
\(w=6/(p-1)\).
[Voight, Example 15.5.7 and Theorem 25.1.1](sources/voight-quaternion-algebras-author.pdf)

Use the all-integral-ideal Hecke operators on the space of all
right-ideal classes, with auxiliary level when needed. The coarse
constant eigenvalue is \(\sigma(m)\) for \(p\nmid m\).
On a one-prime fine quotient the constant and quadratic norm
characters have eigenvalues \(\sigma(m)\) and
\((m/q)\sigma(m)\), with \(q\) the auxiliary prime, and \(m\)
coprime to the auxiliary prime and \(p\). At \(q=p\) use the
radical quotient and require \(p\nmid m\).
Hasse–Schilling gives the positive rational norm image. The
reduction kernel restricts the local norms to square units at the
auxiliary prime, while maximal units elsewhere have all unit norms.
This gives precisely the stated norm characters. Strong approximation
excludes further constituents with one-dimensional factors at a
split place. Jacquet–Langlands transfers the remaining constituents
at trivial definite infinite type to weight two. At each Hecke
index prime the local level is maximal and the algebra is split.
The all-ideal recurrence and the weight-two Ramanujan bound give
operator norm at most \(\tau(m)\sqrt m\) on the complement.
These normalizations are the same for every \(p\).
[All-class normalization and its local hypotheses](PRIMITIVE_NORM_HEEGNER_CORE_SUMSETS.md#quaternionic-sieve) ·
[Classical Hecke transfer](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex)

For \(n\in L_p\), keep the raw shell when 2 splits, with
\[
 S=\sigma(n),\qquad \Delta=\tau(n)\sqrt n .
\]
When 2 is inert, write \(n=2^t m\) with \(t\ge1\) and \(m\) odd.
In the even shell the two block norms are both odd or both blocks
have even coordinates, by anisotropy of the field norm at 2.
The filter \(T_n-T_{n/4}\) keeps the former shell, with the second
term zero when its index is nonintegral. Its constant and
complementary bounds are
\[
 S=3\cdot2^{t-1}\sigma(m),\qquad
 \Delta=2^{t/2}(t+3)\tau(m)\sqrt m/2 .
\]
The dyadic Satake sum has two endpoint coefficients 1 and its
\(t-1\) interior coefficients become \(1/2\) after subtraction.
Scalar multiplication by 2 is projectively trivial at the odd
auxiliary primes, so the filter also counts the marked odd shell.
For both cases,
\[
             S\ge n,\qquad E:=\Delta/S\le\tau(n)/\sqrt n .      \tag{2}
\]
In the dyadic case use
\(E=(t+3)\tau(m)m/(3\sqrt n\,\sigma(m))\) and
\((t+3)/3\le t+1\).

The coarse principal-return count has lower bound
\(4[wS-(1-w)\Delta]\). The factor 4 is the number of norm-one
quaternion units. The correlation of the principal class indicator
with itself has constant mass \(w\), and its complementary contribution
at that class has absolute value at most \((1-w)\Delta\).

For an odd inert \(q\nmid n\), a bad point has one zero block modulo
\(q\). Its principal-fiber image is the normalizer of the nonsplit
torus in \(\mathrm{PGL}_2(\mathbb F_q)\), with relative mass
\(\mu_q=2/(q(q-1))\). Its indicator has mean \(w\mu_q\) on the
fine quotient. The quadratic determinant character averages to
zero on the normalizer, because the torus norm onto
\(\mathbb F_q^\times\) is surjective. Right averaging over the
normalizer turns evaluation at the base point into the normalized
self-correlation. Thus the bad count divided by four is at most
\(w\mu_q S+\Delta\). The measure is the full inverse-unit
measure, and every fine stabilizer is trivial. A noncentral unit has
trace discriminant \(-3\) or \(-4\), while scalar reduction of its
matrix at a split place \(q\) forces divisibility by \(q^2\).
This includes \(q=3\).

For an odd inert \(q\mid n\), one zero block forces both blocks
zero modulo \(q\). At valuation one there are no bad points.
At valuation at least two the bad count is the complete filtered
shell at \(n/q^2\). Here \(S(n/q^2)\le S(n)/q^2\) and
\(\Delta(n/q^2)\le \Delta(n)/q\), by the divisor formulas in
(2). The same bound \(w\mu_q S+\Delta\) is therefore valid without
invoking a Hecke operator at its level prime.

At \(p\) use the radical quotient
\(\mathbb F_{p^2}^\times/\mathbb F_p^\times\). Scalar reduction of a
noncentral unit would force \(p\mid-3\) or \(p\mid-4\), impossible.
The two bad axes, where one block norm vanishes modulo \(p\),
give a single point in the principal fiber after dividing by the
projective units. Its mean is \(\nu=2w/(p+1)\) on the fine space.
Its projections onto the constant and quadratic norm characters
have coefficient \(\nu\) each. The bad count divided by four is
at most \((1+(n/p))\nu S+(1-2\nu)\Delta\) from the same
self-correlation. In particular the main loss is at most
\(4wS/(p+1)\). Removing the axes also removes all zero blocks.
[Local normalizer and radical correlations](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md#small-bad-prime-counts-and-the-ramified-filter)

### The large-prime bound without class number one

The \(Q_p\)-vectors of positive norm \(v\) give principal integral
ideals of norm \(v\), with at most two vectors for each ideal.
The number of all integral ideals of norm \(v\) is at most
\(\tau(v)\), by ideal factorization at inert, split and ramified
primes. Thus \(r_{Q_p}(v)\le2\tau(v)\). This is the only binary
representation bound needed. For example, at \(p=23\) the prime
2 splits, while \(x^2+xy+6y^2\) cannot have value 2. The counting
argument still applies to its principal norm shell. A surviving positive
principal norm belongs to
\(C_p\) just because both inert and ramified primes were excluded.

Let \(H_p(X)\) count nonzero binary vectors with norm at most \(X\).
From \(4Q_p=(2x+y)^2+py^2\), the row at \(y\) has length
\(2\sqrt{X-py^2/4}\) when \(|y|\le2\sqrt{X/p}\).
Its integer count is at most the length plus one. Monotonicity
on the half-line changes the sum of row lengths from the ellipse
area by at most \(2\sqrt X\), and the row count after removing the
origin contributes at most \(4\sqrt X/\sqrt p\). Consequently
\[
 H_p(X)\le \frac{2\pi X}{\sqrt p}+
                  (2+4/\sqrt p)\sqrt X
          \le \frac{7X}{\sqrt p}+4\sqrt X.                    \tag{3}
\]
The factor \(p^{-1/2}\) in the area term is used below.

For an odd inert \(q>P\), count bad points surviving the ramified
filter, so each block is nonzero. A block with norm divisible by \(q\)
is \(q\) times a nonzero binary vector of norm at most \(n/q^2\).
The other block has at most \(2\max_{v\le n}\tau(v)\) representations.
Counting both orientations, the bad count divided by four is
bounded by \(H_p(n/q^2)\max_{v\le n}\tau(v)\). Such a bad prime
cannot exceed \(\sqrt n\).

Summing \(2/(q(q-1))\) over all odd integers \(q\ge3\) gives
\(2(1-\log2)<5/8\), since the first two terms of the
\(\operatorname{arctanh}(1/3)\) series give
\(\log2>56/81>11/16\).
Use \(\#\{q\le P\}\le P\),
\(\sum_{q>P} q^{-2}\le 2/P\) for \(P\ge1\), and
\(\sum_{q\le\sqrt n}q^{-1}\le1+\tfrac12\log n\). These sums may
all be enlarged to integer sums in the stated ranges.
The total loss beyond the two main losses \(4/(p+1)\) and \(5/8\),
as a fraction of \(4wS\), is at most
\[
 \frac{\max_{v\le n}\tau(v)}{w}
 \left(\frac{P}{\sqrt n}+\frac{14}{\sqrt p\,P}
                           +\frac{6+2\log n}{\sqrt n}\right).   \tag{4}
\]
Here (2) supplies the spectral term, and (3) supplies the two
large-prime terms. The first two main losses leave more than
\(1/24\) for \(p\ge11\).

Take \(P=(n/p)^{1/4}\). For \(n\ge p^2\) we have
\(\log n\le8n^{1/8}\) and \(p^{1/4}\le n^{1/8}\).
Using \(1/w<p/6\), (4) is then bounded above by
\[
                    \frac{37}{6}\max_{v\le n}\tau(v)
                                    \frac{p^{3/4}}{n^{1/4}}.   \tag{5}
\]
The exponent \(3/4\) comes from the balance between the
spectral term and the discriminant-sensitive lattice area term.

### Explicit divisor bounds

For \(q\ge2^h\), the Euler factor \((a+1)^h/q^a\le1\) because
\(a+1\le2^a\). At smaller primes the ratio of successive factors
for \(a\ge2h\) is at most
\((1+1/(2h+1))^h/2<1\).
For the last inequality bound the binomial sum by the geometric
series with common ratio \(h/(2h+1)\). Thus the finite product
defining \(T_h\) gives \(\tau(v)\le T_h v^{1/h}\) for every \(v\).
At (1) the bound in (5) is at most \(37/(6\cdot256)<1/24\).
This leaves a positive shell and proves the theorem.

The exact rational certificate gives the following three
choices. They describe thresholds for all eligible targets,
including the earlier cases \(p=3,7\).

| \(h\) | Upper bound for \(T_h\) | Explicit upper bound for \(M(p)\) |
| --- | ---: | ---: |
| 8 | \(2^{16}\) | \(2^{192}p^6\) |
| 10 | \(2^{36}\) | \(2^{294}p^5\) |
| 16 | \(2^{685}\) | \(2^{3696}p^4\) |

For any fixed \(\varepsilon>0\), choose an integer
\(h\ge4+12/\varepsilon\) to get the effective exponent
\(3h/(h-4)\le3+\varepsilon\).
The finite splitting law uses prime equidistribution only with
the window fixed. The polynomial existence proof uses no
prime-distribution estimate in the variable discriminant.

The preceding [fixed-modulus asymptotic](../../results/pointwise-binary-splitting-semigroups/theorems.tex)
gives gap-free tails for fixed splitting data in its stated scope.
Here the quaternionic proof supplies a different kind of uniformity,
an unconditional existence threshold polynomial in \(p\).
The [five exact cases](PRIMITIVE_NORM_HEEGNER_CORE_SUMSETS.md)
supply much smaller thresholds at the benchmark discriminants.

## Reproduction

From the repository root, run with Python assertions active.

```sh
python3 research/arithmetic/verify_prime_discriminant_core_uniform.py
```

The checker uses integer prime sieving and rational comparisons
for the three uniform divisor products. It also recovers the
complete Bernoulli pattern distribution through \(B=12\) and
checks equal reduced-residue counts modulo \(8\cdot3\cdot5\cdot7\cdot11\).
The retained receipt records exact values and input hashes.
[Certificate and CRT control](prime_discriminant_core_uniform_checks.json)
