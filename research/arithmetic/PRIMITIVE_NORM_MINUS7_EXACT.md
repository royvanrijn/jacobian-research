# Exact sums of two primitive discriminant -7 norms

Let
\[
 Q(x,y)=x^2+xy+2y^2,\qquad
 P=\{Q(x,y)>0:\gcd(x,y)=1\}.
\]
Signs and zero coordinates are allowed. Let \(C\) contain 1 and the
positive integers supported on primes congruent to \(1,2,4\pmod7\);
in particular 2 is allowed.

**Theorem M7.** The exact identities are
\[
 C+C=\{n\ge2:7\nmid n\}\setminus\{11,29\},                 \tag{M7.1}
\]
\[
 P+P=\{N\ge2:49\nmid N\}\setminus\{7,77,203\}.             \tag{M7.2}
\]
Moreover, for every \(N\ge2^{69}\) with \(49\nmid N\), the number
of ordered pairs of positive primitive norm values satisfies
\[
       \#\{(a,b)\in P^2:a+b=N\}\ge N^{2/3}/20000.          \tag{M7.3}
\]

The proof uses a quaternion norm shell, a simultaneous ramified/inert
prime filter, exact uniform interval certificates, and a descent that
handles exceptional prime quotients. The finite scan is not used as
evidence for an unbounded extrapolation. Its upper endpoint is joined
to a proved analytic range. This is an internally reviewed deduction
from classical inputs; no formal verification, external review or
historical priority is claimed.

## Primitive norms and ramification

Put \(R=\mathbb Z[\omega]\), \(\omega=(1+\sqrt{-7})/2\), and
\(\pi=\sqrt{-7}\). The Minkowski bound
\((2/\pi)\sqrt7<2\) shows that \(R\) has ideal class number one:
every ideal class has an integral representative of norm one.
Primitive ideal factorization therefore gives the disjoint union
\[
                           P=C\ \dot\cup\ 7C.             \tag{M7.4}
\]
An inert prime dividing a norm divides both coordinates. Ramified
valuation is at most one for a primitive element. Conversely, orienting
each entire split prime power to one of its principal prime ideals,
and optionally multiplying by \(\pi\), gives every value in (M7.4).
This proves primitive existence for a value; the coordinate witnesses
used later need not themselves be primitive at split primes.

The matrix of multiplication by \(\pi\) on the basis \(1,\omega\) is
\[
 M=\begin{pmatrix}-1&-4\\2&1\end{pmatrix},\quad
 M^2=-7I,\quad\det M=7,\quad Q(Mv)=7Q(v).
\]
It maps 7-unit primitive generators to primitive generators of norm
seven times the original norm. At 7 the output can have both
coordinates divisible only if \(2x+y\equiv0\pmod7\), which is exactly
the condition that the original norm is divisible by 7.
Thus this is not an assertion that all of \(P\) is closed under
multiplication by 7.

Modulo 7, \(Q(x,y)=(x+4y)^2\). Primitive 7-unit norms are
nonzero squares, and primitive ramified norms have valuation one
and square unit quotient. Since \(-1\) is nonsquare, a sum
divisible by 7 must use two ramified summands, and a sum cannot
have valuation at least two. In particular
\[
 P\cap7\mathbb Z=7C,\qquad
             7m\in P+P\ \Longleftrightarrow\ m\in C+C.     \tag{M7.5}
\]
These identities hold for every positive \(m\). Multiplication of
two \(C\)-summands by any split prime preserves their membership
in \(C\).

## The maximal quaternion order and spectral input

Let
\[
 B=(-7,-1)_{\mathbb Q},\qquad O=R\oplus Rj,\quad
 j^2=-1,\quad ja=\bar a j.
\]
Its norm is \(Q\oplus Q\). Voight, Theorem 15.5.5 and
Example 15.5.7, give algebra discriminant 7 and maximality of this
order. Its units are \(1,-1,j,-j\): the only field elements of
norm one are \(1,-1\), and the two nonnegative field norms sum
to one. The projective unit group \(\Gamma\) has order two.

Voight's maximal-order mass formula, Theorem 25.1.1, gives total
projective mass \((7-1)/12=1/2\). The principal class already
contributes \(1/|\Gamma|=1/2\); hence the entire ideal-class set
has one element. All integral right ideals of norm \(n\), for
\(7\nmid n\), are principal and have four generators. The all-ideal
degree formula consequently gives the full signed ordered shell
count
\[
 r_O(n)=\#\{\alpha\in O:\operatorname{nrd}\alpha=n\}
                              =4\sigma(n)\quad(7\nmid n). \tag{M7.6}
\]
Here \(\sigma\) is the divisor sum. The same ideal factorization
in \(R\) gives, for \(t\ge1\),
\[
 r_Q(t)=2\sum_{d\mid t}\chi_{-7}(d),\qquad
 r_Q(t)\le2M(X)\quad(1\le t\le X),
\]
where
\[
 M(X)=\max_{1\le t\le X}
             \prod_{\substack{\ell\mid t\\\chi_{-7}(\ell)=1}}
                           (v_\ell(t)+1).                \tag{M7.7}
\]

We specify the spectral input and its application. At an odd inert
prime \(q\), reduce to \(G_q=\mathrm{PGL}_2(\mathbb F_q)\).
At 7 use the radical quotient
\(G_7=\mathbb F_{49}^{\,*}/\mathbb F_7^{\,*}\), cyclic of order eight.
For one or both of these auxiliary levels the complete quotient is
\(\Gamma\backslash\prod G_\ell\), with uniform probability measure.
The unit image injects: at odd \(q\), a scalar matrix of trace zero
cannot also have norm one; at 7, \(j\notin\mathbb F_7\).
This includes \(q=3\).

For an index \(n\) coprime to the auxiliary levels and 7, the all-ideal
operator is self-adjoint, commutes with the right finite group action,
and at the base point satisfies
\[
 (T_nh)(x_0)=\frac14
       \sum_{\operatorname{nrd}\alpha=n}h(\rho(\alpha)).
                                                               \tag{M7.8}
\]
Principalization initially gives inverse generators. Projectively
\(\alpha^{-1}=\bar\alpha\), and conjugation permutes the complete
shell, giving the displayed orientation and factor.
Voight, Propositions 41.3.1 and 41.3.6, specify the degree
\(\sigma(n)\) and the all-integral-ideal prime-power recurrence.

At level 7 the global norm characters are \(1,\chi_7\); at joint
level \(q,7\) they are
\[
                        1,\chi_7,\chi_q,\chi_7\chi_q .
\]
Here \(\chi_q\) is the determinant square character and \(\chi_7\)
the residue-field norm square character. These are the complete
one-dimensional **global automorphic** constituents. Reduced norms of maximal units give all
units away from the auxiliary levels, while the reduction kernels
and central scalars leave only the two local square-class choices.
Rational reduced norms are all positive rationals by the
Hasse--Schilling theorem. Strong approximation for the norm-one
group at a split place excludes additional constituents that are
characters at such a place. This is the explicit single- and
two-level extension of the argument in the
[retained spectral proof](../../results/eisenstein-restricted-prime-sums/sections/prime-bounds.tex).

For precision, this exclusion concerns the full smooth local factor
\(\pi_\ell\) at a split finite place, rather than its finite-dimensional
space of level-fixed vectors. If \(\pi_\ell\) is finite-dimensional,
the intersection of the stabilizers of a basis is an open normal
kernel. Conjugating small upper and lower unipotents by diagonal
matrices then puts all of \(\mathrm{SL}_2(\mathbb Q_\ell)\) in the
kernel. Given a vector \(f\), \(x\in B^\times(\mathbb A_f)\), and
\(h\in B^1(\mathbb A_f)\), norm-one strong approximation with
\(S=\{\infty,\ell\}\) supplies \(\gamma\in B^1(\mathbb Q)\) such that
\(k=x^{-1}\gamma^{-1}xh\) lies in the right stabilizer of \(f\) away
from \(\ell\). At \(\ell\), \(k\) has norm one and acts trivially.
Thus \(f(xh)=f(\gamma xk)=f(x)\). The trivial infinite type supplies
archimedean invariance, so the constituent factors through reduced
norm and is one of the global characters just listed. This verifies
the split-place hypothesis of Jacquet--Langlands Theorem 14.4 for
every remaining constituent.

Characters of the finite right action need not be global norm
characters. Already \(\Gamma\backslash G_7\simeq C_4\) has two
quartic characters. The norm-two shell has eight pure-component
points, all reducing to \(\Gamma\), and four mixed points
\(\pm1\pm j\), all reducing to \(h=[1+j]\). Since
\((1+j)^2=2j\) and \(h\ne1\), \(h\) has order two, and (M7.8) gives
\[
                   T_2=2I+R_h,\qquad
                   \operatorname{Spec}(T_2)=\{3,3,1,1\}.
\]
Both global norm characters take \(h\) to \(1\), since
\((2/7)=1\). The quartic finite characters take \(h\) to \(-1\)
and have eigenvalue \(1\); they remain in the Ramanujan complement.

On the complement of those norm characters,
\[
                            \|T_n\|\le\tau(n)\sqrt n.      \tag{M7.9}
\]
The imported argument is global Jacquet--Langlands to weight two,
followed by the weight-two Ramanujan bound. Every prime dividing
the index is split in \(B\) and has maximal local level. Its
all-ideal eigenvalue at \(\ell^a\) is a sum of \(a+1\) monomials
in Hecke roots of absolute value \(\sqrt\ell\); multiplication
gives (M7.9). The source ledger pins the actual Voight, Jacquet--Langlands
and Deligne editions and the existing proof bytes. No operator
bound at a prime dividing its added level is assumed.
The prime 2 is split in \(R\) and remains at maximal quaternion
level, so arbitrary powers of 2 cause no new level condition.

## A simultaneous filter with the correct normalization

For \(7\nmid n\), put
\[
 E(n)=\frac{\tau(n)\sqrt n}{\sigma(n)},\qquad
              \gamma(n)=\frac{3-\chi_{-7}(n)}4.
\]
Let \(G(n)\) count shell points whose two field norms are 7-units.
Both norms are then positive. On the four-point quotient
\(\Gamma\backslash G_7\), the two component axes are the single
base point: their projective residues are \(1,j\).
Its indicator \(\delta\) has mean \(1/4\) and the same coefficient
on \(\chi_7\); its complementary squared norm is \(1/8\).
Evaluating its correlation and using (M7.9) gives
\[
 \frac{G(n)}{4\sigma(n)}\ge\gamma(n)-\tfrac12E(n).         \tag{M7.10}
\]
If \(\chi_{-7}(n)=-1\), no axis can occur at all, since the other
norm would be a nonzero square modulo 7. In that case
\(G(n)=4\sigma(n)\) exactly.

For an inert \(q\nmid n\), let \(H_q\) be the nonsplit torus
normalizer. It consists precisely of the invertible residues with
one field component zero. It contains \(\Gamma\), has order
\(2(q+1)\), and fraction
\[
                         \mu_q=\frac{2}{q(q-1)}
\]
in \(G_q\). Its determinant square character has mean zero, since
the torus norm is onto \(\mathbb F_q^*\).

On the joint quotient define two indicators:
\(u\) requires the \(q\)-coordinate in \(H_q\) and the 7-coordinate
in \(\Gamma_7=\{1,j\}\); \(v\) requires the \(q\)-coordinate in
\(H_q\) and the 7-coordinate outside \(\Gamma_7\).
They are right-invariant under \(H_q\times\Gamma_7\); the base
orbit is precisely the support of \(u\). Therefore
\[
 (T_nv)(x_0)=\frac{\langle u,T_nv\rangle}{\mathbb E u}.
\]
This identity does not require \(v(x_0)=1\).

Write \(a=\mu_q/4\). The means of \(u,v\) are \(a,3a\);
their \(\chi_7\) coefficients are \(a,-a\), and both other
character coefficients vanish. Their complementary squared norms
are \(a-2a^2\) and \(3a-10a^2\).
If \(b_q(n)\) counts points bad at \(q\) inside the 7-unit shell,
Cauchy's inequality gives
\[
 \frac{b_q(n)}{4\sigma(n)}
 \le\gamma(n)\mu_q+
 \sqrt{(1-\mu_q/2)(3-5\mu_q/2)}\,E(n)
 \le\gamma(n)\mu_q+\tfrac74E(n).                         \tag{M7.11}
\]
Thus the ramified filter is retained in the main term of the
bad-prime count.

If \(q\mid n\), inert anisotropy forces both coordinate blocks
divisible by \(q\) as soon as either is. The count is zero at
valuation one; otherwise it is exactly \(G(n/q^2)\). In particular
\[
 \frac{b_q(n)}{4\sigma(n)}
 \le\frac{\sigma(n/q^2)}{\sigma(n)}<q^{-2}
                               \le\gamma(n)\mu_q.        \tag{M7.12}
\]
This handles every valuation without applying a level-\(q\) operator
to an index divisible by \(q\). In the negative-character case
one may omit level 7 entirely. The ordinary \(H_q\) correlation
then has complementary norm \(\mu_q(1-\mu_q)\), giving
\(b_q/(4\sigma)\le\mu_q+(1-\mu_q)E\le\mu_q+E\).

## Large bad primes

For every real \(X\ge0\), the number of nonzero \(R\)-vectors with
norm at most \(X\) is at most \(3X\). Indeed, row integration gives
\[
 S_Q(X)\le\frac{2\pi}{\sqrt7}X+
                           (2+4/\sqrt7)\sqrt X.
\]
For integer \(X\ge36\), use \(\pi<22/7,\sqrt7>66/25\) to obtain
\(S_Q(X)/X\le50/21+58/99<3\).
The checker completely enumerates \(1\le X\le35\) and verifies
the same bound; flooring extends it to real \(X\).

For an inert \(q\) dividing one component norm, that component is
\(q\) times a nonzero \(R\)-vector. Count it by \(3n/q^2\), the
other component by \(2M(n)\), and both orientations. Hence
\[
                 b_q(n)\le12M(n)n/q^2.                  \tag{M7.13}
\]
Only \(q\le\sqrt n\) can occur.

Let \(\mathcal I\) be the inert primes and define
\[
 B=\sum_{q\in\mathcal I}\mu_q,\quad
 \Pi(z)=\#\{q\in\mathcal I:q\le z\},\quad
 T(z)=\sum_{\substack{q\in\mathcal I\\q>z}}q^{-2}.
\]
Combining (M7.10)--(M7.13), the surviving shell fraction is at least
\[
 \gamma(n)(1-B)-E(n)\{\tfrac12+\tfrac74\Pi(z)\}
                                      -3M(n)T(z).        \tag{M7.14}
\]
Strict positivity gives \(n\in C+C\): surviving positive norms
have no inert factors and are 7-units, and (M7.4) supplies
primitive representations of their values. In the branch
\(\chi_{-7}(n)=-1\) the sharper sufficient inequality is
\[
                  E(n)\Pi(z)+3M(n)T(z)<1-B.              \tag{M7.15}
\]

## Exact arithmetic envelopes

The inert odd primes lie in \(3,5,13\pmod{14}\).
Counting these three progressions over a complete period and using a
first term plus a decreasing integral gives, for \(z>0\),
\[
 \Pi(z)\le3z/14+13/14,\qquad
                   T(z)\le3/(14z)+3/z^2.                 \tag{M7.16}
\]
The constant \(13/14\) is the maximum discrepancy of the three
progressions within one period. The weaker constant \(3/2\) is
also useful below.

For a retained prime-table limit \(R\), sum the actual inert
\(\mu_q\) through \(R\), rounding upward at scale \(2^{64}\).
For the rest use
\(\mu_q\le2(R+1)q^{-2}/R\) and the three progressions starting
just above \(R\). This gives a rational bound \(\bar B\) for
\(B\), and the certificates prove
\[
                              B\le\bar B<47/100.          \tag{M7.17}
\]
No floating point comparison enters the proof.

The complete Euler-product certificates give
\[
 M(X)\le8X^{1/6}\quad(X\ge1),\qquad
 E(n)<\frac{91}{20}n^{-1/3}
           \quad(n\text{ supported on inert primes}).    \tag{M7.18}
\]
Here is the completeness argument. For \(M\), each normalized
sixth-power factor is \((a+1)^6/\ell^a\).
For \(\ell\ge64\) it is at most one because \(a+1\le2^a\).
At smaller split primes, exponents \(0,\ldots,8\) suffice:
after \(a=8\) the next ratio is at most \((10/9)^6/2<1\).
Multiplication of the exact finite maxima gives a product less
than \(8^6\).

For \(E(n)n^{1/3}\), the sixth-power factor is
\[
                 W_\ell(a)=
     \frac{(a+1)^6\ell^{5a}}{\sigma(\ell^a)^6}.            \tag{M7.19}
\]
Its consecutive ratio is strictly less than
\(((a+2)/(a+1))^6/\ell\).
At the inert primes below 64, it is therefore enough to use
\(a\le4\) at 3, \(a\le3\) at 5, and \(a\le1\) at
\(13,17,19,31,41,47,59,61\).
At all larger primes the factor is at most one.
The exact product of the finite maxima is less than \((91/20)^6\).

For bounded endpoints the checker uses stronger envelopes.
A maximizing \(M\)-value is obtained by placing nonincreasing
exponents on an initial segment of the smallest split primes.
Exchanging out-of-order exponents reduces the integer cost and
preserves its divisor product. A complete finite recursion through
\(2^{54}\) has 2556 states; the recorded prime-prefix product and
exponent bound prove completeness.

For \(E\), enumerate all 5120 tuples within the caps in (M7.19).
For a bound \(U\), retain the tuples of integer cost at most \(U\),
then maximize their sixth-power products. Truncating any larger
actual exponent to its cap increases the factor and reduces its
cost. Removing primes at least 64 does the same. This proves that
the finite frontier bounds every inert-supported target up to \(U\).
Split the tuples into the nine categories
\[
             (v_3(n),v_5(n))\in\{0,1,\ge2\}^2.
\]
The truncation preserves the category. Upward sixth roots are
taken at denominator \(10^6\) and checked by sixth powers.
If the resulting constant is \(c_{a,b}(U)\), then on \(L\le n\le U\),
\[
                  E(n)\le
       c_{a,b}(U)/\lfloor\sqrt[3]L\rfloor.                \tag{M7.20}
\]

## Finite point counts close the intermediate range

When \(\chi_{-7}(n)=1\), the ordered pair of nonzero square
residues summing to \(n\pmod7\) is unique: both are \(n/2\).
If the first block is \(q\alpha\), then
\(Q(\alpha)\equiv n/(2q^2)\pmod7\), one nonzero square class.
For \(r=1,2,4\), let
\[
 S_r(X)=\#\{(x,y)\in\mathbb Z^2:0<Q(x,y)\le X,\
                                      Q(x,y)\equiv r\pmod7\},
 \qquad F(X)=\max_{r=1,2,4}S_r(X).
\]
These are ordinary lattice counts, without a gcd restriction.

The checker enumerates the complete table through \(X_0=2^{18}\).
For \(y\ge0\), put \(h=\lfloor\sqrt{4X_0-7y^2}\rfloor\);
every possible \(x\) lies in
\[
 -\lfloor(h+y)/2\rfloor\le x\le\lfloor(h-y)/2\rfloor,
 \quad 0\le y\le\lfloor\sqrt{4X_0/7}\rfloor.
\]
Weight a row by two for \(y>0\), to include its simultaneous
sign-reversed row, and by one at \(y=0\).
This is a complete enumeration of the signed ordinary coefficients.
Cumulative sums give every \(S_r\) and \(F\). A separate
divisor-character computation checks the coefficients.
The first seven \(F\)-values are \(2,4,4,6,6,6,6\).

For a fixed interval \([L,U]\), choose an initial segment of the
inert primes as the small-prime set. The checker verifies that
every remaining inert \(q\le\sqrt U\) satisfies
\(\lfloor U/q^2\rfloor\le X_0\).
It includes 3 and 5 in the small-prime set. If \(j\) is its size,
the large-prime fraction is bounded exactly by
\[
 \frac{M(U)\rho_{a,b}}{L}
       \sum_{\substack{q\in\mathcal I\text{ outside the cut}\\q\le\sqrt U}}
                         F(\lfloor U/q^2\rfloor),        \tag{M7.21}
\]
where
\[
 \rho_{a,b}=
 \frac{3^a}{\sigma(3^a)}\frac{5^b}{\sigma(5^b)},\qquad
                         a,b\in\{0,1,2\}.
\]
The representative exponent 2 gives an upper bound for every
larger exponent, since \(\ell^a/\sigma(\ell^a)\) decreases.
Thus \(n/\sigma(n)\le\rho_{a,b}\).
The factor in (M7.21) is \(M\), not \(2M\):
two choices of bad block times the coefficient bound \(2M\),
divided by \(4\sigma(n)\), give \(M/\sigma(n)\).

For a known divisor \(p=3\) or 5 of \(n\), omit its small-prime
spectral error and refund \(\gamma\mu_p-\delta_p\) to the baseline,
where \(\delta_p=0\) at valuation one and \(\delta_p=p^{-2}\)
at larger valuation. This follows from (M7.12); removing the
\(\mu_p\) contribution without retaining \(\delta_p\) would be invalid.
Let \(d\) be the number of these known divisors. In the positive
character case the finite sufficient inequality is
\[
 e\{\tfrac12+\tfrac74(j-d)\}+\text{large tail}
 <\tfrac12(1-\bar B)
                  +\sum_{\substack{p=3,5\\p\mid n}}(\mu_p/2-\delta_p).
                                                               \tag{M7.22}
\]
Here \(e\) is (M7.20). In the negative character case the small
error is \(e(j-d)\), the baseline is
\(1-\bar B+\sum(\mu_p-\delta_p)\), and the ordinary bound
(M7.13) gives the large term
\[
 3M(U)\rho_{a,b}
       \sum_{\substack{q\text{ outside the cut}\\q\le\sqrt U}}q^{-2}.
\]
Reciprocal squares are rounded upward at scale \(2^{64}\).

The complete certificate partitions \([2^{24},2^{37})\) into
257 adjacent integer intervals. It uses sixteen equal subdivisions
per dyadic interval and further cuts at envelope jumps.
On every interval it checks all eighteen combinations of character
sign and the nine valuation categories. It also checks one
positive-character boundary family using
\(E(n)\le4/\lfloor\sqrt L\rfloor\), valid whenever \(\tau(n)\le4\).
There are exactly 4883 cases. Every case has strict margin greater
than \(1/64\); no uncovered valuation or interval remains.
The worst interval is
\[
 L=17825792,\quad U=18874367,\quad
 \chi_{-7}(n)=1,\quad v_3(n)=v_5(n)=0,
\]
with \(M(U)=112\), \(e=1907201/261000000\), and cut at 73.
Its total error is
\[
                 \frac{563281558829}{2271744000000},
\]
below its certified baseline by more than \(1/64\).
All primes needed for these finite tails are enumerated through
\(2^{19}>\sqrt{2^{37}}\). Every point-table query is explicitly
checked within its certified range. Recording the last included
prime rather than a real cutoff causes no issue: it is the
included prime set and the verified query condition that matter.

## The unbounded tail and the prime boundary

For \(k=37,\ldots,53\), exact divisor and \(E\)-frontiers supply
\(d_k=M(2^{k+1}-1)\), \(c_k\ge1\) and
\(e_k=c_k/\lfloor\sqrt[3]{2^k}\rfloor\).
The checker verifies all seventeen inequalities
\[
 e_k\{\tfrac12+\tfrac74\Pi(z_k)\}
                    +3d_kT(z_k)<(1-\bar B)/2-1/64.       \tag{M7.23}
\]
Actual inert primes and the progression tail in (M7.16) give the
rational upper bound for \(T(z_k)\). These inequalities cover
every inert-supported target in the stated intervals.

They cover the prime boundary as well. For \(\tau(n)\le4\),
\(E(n)\le4/\sqrt n\). At these seventeen lower endpoints the
checker explicitly verifies
\[
      4\lfloor\sqrt[3]{2^k}\rfloor\le\lfloor\sqrt{2^k}\rfloor.
\]
Thus \(4/\sqrt n\le e_k\) throughout each interval.
This is a separate domination argument; no inert-support Euler
bound is silently applied to split primes.

For inert-supported \(n\ge2^{54}\), use (M7.18) and take
\(z=(7/4)n^{1/4}\). The error in (M7.14) is at most
\[
 \frac{185799}{31360}n^{-1/12}
           +\left\{\frac{17}{8}\frac{91}{20}
                                      +\frac{1152}{49}\right\}n^{-1/3}.
\]
Since \((45/2)^{12}<2^{54}\) and \(n^{1/3}\ge2^{18}\),
this is below \(53/200\), by the exact margin
\[
                         \frac{9574319}{6165626880}>1/1024.
\]
The baseline is greater than \(53/200\) by (M7.17).
For every \(\tau(n)\le4\) target,
\(E(n)\le4n^{-1/2}\le(91/20)n^{-1/3}\), so the same
unbounded tail applies to it. Consequently every inert-supported
target, and every split-prime target, above \(2^{24}\) is in \(C+C\).

## Finite coverage and complete exceptional descent

The checker constructs the full primitive value flag through
\(B_0=2^{24}\) in two ways: complete gcd-restricted coordinates,
and the character-factor criterion (M7.4).
They agree at every index. Deleting multiples of 7 gives the
complete \(C\)-flag. There are 2490663 primitive values and
2161776 \(C\)-values.

Packed shifts by all positive primitive anchors at most 5000
cover every \(49\)-admissible target through \(B_0\) except
\(7,77,203\). The 1050 anchors give only genuine sums.
Likewise the 903 \(C\)-anchors cover every 7-unit target in
\([2,B_0]\) except \(11,29\).
Every remaining target is checked by a complete positive-pair
search through half the target; there are no unresolved positive
residuals. All local exclusions are also checked, and the two
ramified masks agree with (M7.5).
Positivity ensures no value above \(B_0\) can contribute to a
target in this finite range.

Suppose now that \(n\) is a least exception to (M7.1) outside
\(11,29\). It exceeds \(B_0\). If a split prime \(q\mid n\),
minimality and multiplicative closure force
\(n/q\in\{1,11,29\}\).
The first case makes \(n\) a split prime, already covered by the
analytic boundary argument. In the other cases \(n=11q\) or
\(29q\); dividing by the other split factor forces
\(q\in\{11,29\}\), or else a smaller represented target lifts.
The only remaining products are \(121,319,841\), already covered;
explicit witnesses are
\[
                121=8+113,\quad319=2+317,\quad841=32+809.
\]
Every displayed summand belongs to \(C\).
If \(n\) has no split factor, all its prime factors are inert,
and the analytic core bounds give the final contradiction.
This proves (M7.1). Conversely, no two nonzero square residues
sum to zero modulo 7, and the finite negative checks prove
the two stated unit exceptions. For instance
\[
 C\cap[1,29]=\{1,2,4,8,11,16,22,23,29\}
\]
makes their nonrepresentation direct.

For \(7\nmid N\), (M7.1) gives \(N\in P+P\) unless
\(N=11,29\); these have the primitive patches
\(11=4+7\) and \(29=1+28\).
For \(7\mid N\) with \(49\nmid N\), use (M7.5).
Its quotient fails exactly at \(1,11,29\), including the positive
summand obstruction at quotient 1. This gives exactly
\(7,77,203\). Multiples of 49 are already excluded locally.
The proof of (M7.2) is complete.

## An explicit count at every sufficiently large target

A further complete Euler certificate uses (M7.19) over all primes
below 64 except 7. The caps are 8 at 2, 4 at 3, 3 at 5, 2 at 11,
and 1 at every larger prime below 64. The same ratio argument
proves completeness; the exact product is less than \(15^6\).
Therefore
\[
                    E(n)\le15n^{-1/3}\quad(7\nmid n).    \tag{M7.24}
\]
This bound has no inert-support hypothesis.

Take \(z=n^{1/4}\) in (M7.14), use \(M\le8n^{1/6}\) and
the weaker \(\Pi(z)\le3z/14+3/2\). The error is at most
\[
                 \frac{603}{56}n^{-1/12}
                                      +\frac{951}{8}n^{-1/3}.
\]
For \(n\ge2^{66}\), use \(n^{1/12}>45\), \(n^{1/3}\ge2^{22}\).
The error is at most \(281051653/1174405120\), leaving margin
\(150828519/5872025600>1/40\). There are at least \(n/10\)
surviving coordinate points.

A fixed ordered value pair has at most
\((2M(n))^2\le256n^{1/3}\) coordinate witnesses.
Thus there are at least \(n^{2/3}/2560\) ordered \(C\)-pairs.
For \(N\ge2^{69}\) with \(49\nmid N\), write \(N=7^e n\),
\(e=0,1\); then \(n>2^{66}\). Scaling value pairs injectively
into \(P^2\), and using \(7^{2/3}<4\), proves (M7.3).

## Evidence and remaining scope

The [checker](verify_primitive_norm_minus7.py) deterministically
reconstructs the complete finite value masks and ordinary point array.
The [receipt](primitive_norm_minus7_checks.json) retains their hashes
and counts, the Euler/divisor frontiers, all interval cases and exact
global constants; it does not retain the raw mask or point-array bytes.
The generators, finite domains and rejection controls are explicit.
In particular, membership is checked by exact equality of complete
coordinate and prime-support bytearrays before their hashes are recorded.
The [source ledger](primitive_norm_minus7_source_ledger.json) pins
the classical book and spectral inputs and their application scope.

The infinite spectral theorem is imported and applied in the
written proof; it is not established by the finite checker.
Separate internal reviewers reconstructed the joint quotient,
all character projections, interval completeness, point-array
tails and exceptional descent. This is internal review, not
formal verification or independent external review.

Run from the repository root:

    python3 research/arithmetic/verify_primitive_norm_minus7.py --output research/arithmetic/primitive_norm_minus7_checks.json

Certificate verification requires active Python assertions. The checker rejects optimized execution (-O, -OO, or optimization enabled through PYTHONOPTIMIZE) before generating a receipt.

The exact exception lists are proved globally. The explicit
power lower bound has the stated onset; it is not an unweighted
leading asymptotic. The independent \(C_{12}+C_{12}\) barrier,
other exact discriminants, mixed-field exceptions and conductor
extensions are not resolved by this result.
