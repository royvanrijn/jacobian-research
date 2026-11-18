# Exact sums of two primitive discriminant -11 norms

Let
\[
 Q(x,y)=x^2+xy+3y^2,\qquad
 P=\{Q(x,y)>0:\gcd(x,y)=1\}.
\]
Coordinate signs and zero coordinates are allowed. Values, rather than
coordinate witnesses, are counted in the following statement.

Let \(C\) contain 1 and the positive integers whose prime factors all
satisfy \((p/11)=1\).

**Theorem M11.** The following identities hold without exceptions:
\[
 C+C=\{n\ge2:2\mid n,\ 11\nmid n\},\qquad
 P+P=\{N\ge2:2\mid N,\ 121\nmid N\}.                    \tag{M11.0}
\]
For every even \(N\ge2^{71}\) with \(121\nmid N\), the number of
ordered value pairs also satisfies the explicit bound
\[
 \#\{(a,b)\in P^2:a+b=N\}\ \ge\ \frac{N^{2/3}}{100000}.
 \tag{M11.1}
\]

The proof first reduces a least exception to the inert-core family
\[
 N=11^e n,\quad e\in\{0,1\},\quad
 2\mid n,\quad 2\le n<2^{60},\quad
 p\mid n\Longrightarrow (p/11)=-1.                       \tag{M11.2}
\]
The sharper bounds below cover every such \(n\ge2^{25}\) in 35 complete
dyadic intervals. An exhaustive exact certificate covers every eligible
core target through \(2^{25}\), not only those supported on inert primes.
Together these close (M11.2); it is no longer an open finite obligation.

This is an explicit specialization of the classical quaternionic Hecke
method, with a multi-class mass normalization and a ramified-prime filter.
The general eventual binary-pair theorem is already proved in the
[analytic bridge](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md). Its abstract cusp
bounds alone do not give the numerical constants below. No historical
priority, formal verification or external review is claimed.

## Primitive values and the ramified reduction

Put \(R=\mathbb Z[\omega]\), \(\omega=(1+\sqrt{-11})/2\), and
\(\pi=\sqrt{-11}\). The ideal class number of \(R\) is one: the
imaginary quadratic Minkowski bound is \((2/\pi)\sqrt{11}<3\), whereas
2 is inert, so there is no ideal of norm 2. The only representative
ideal norm required is therefore 1. Let \(C\) contain 1 and all positive
integers supported on primes with \((p/11)=1\). Ideal factorization gives
\[
                         P=C\ \cup\ 11C.                 \tag{M11.3}
\]
For completeness, an inert prime dividing a norm divides both coordinates.
A primitive element has ramified norm valuation at most one, since
\((11)=(\pi)^2\). At a split prime, choosing the entire prime power on
one of the two oriented principal prime ideals produces a primitive
generator. Multiplying such generators and optionally \(\pi\) proves
sufficiency in (M11.3). Thus all primitive values are odd.

Use the definite quaternion algebra and order
\[
 B=(-11,-1)_{\mathbb Q},\qquad
 O=R\oplus Rj,\quad j^2=-1,\quad ja=\bar a j.
\]
Its norm is \(F=Q\oplus Q\). The algebra has discriminant 11 and this
order is maximal: its norm Gram determinant is \(121\), so its reduced
order discriminant is 11. These are also precisely the \(p=11\) cases
of Voight, Example 15.5.7 and Theorem 15.5.5. Its units are
\(\{1,-1,j,-j\}\), because \(Q=1\) has only the two solutions
\((\pm1,0)\).

Modulo 11, \(4Q(x,y)=(2x+y)^2+11y^2\). Since \(-1\) is nonsquare
modulo 11, \(11\mid F(a,b)\) forces both \(a,b\in\pi R\).
Multiplication of both components by \(\pi\) is therefore a bijection
of the complete shells of norms \(n\) and \(11n\). In coordinates,
\[
 (x,y)\longmapsto(-x-6y,2x+y),\quad
 M^2=-11I,\quad \det M=11,\quad Q(Mv)=11Q(v).
\]
Consequently \(121\mid N\) prevents separate coordinate primitivity.
For an eligible target write \(N=11^e n\), \(e=0,1\), \(11\nmid n\).
It suffices to find a shell point of norm \(n\) for which both field
norms are units at 11 and contain no inert prime. Its positive values
are in \(C\); multiplying them by \(11^e\) gives the required members
of \(P\). This stronger unit-unit restriction is possible for every
unit target modulo 11: there are \(10-2(n/11)>0\) ordered nonzero
solutions of \(u^2+v^2=n\). The matrix \(M\) is the identity modulo 2.

## The exact shell and its odd-component filter

Write \(\sigma(t)=0\) if \(t\) is not a positive integer, and let
\[
 f_{11}(\tau)=\eta(\tau)^2\eta(11\tau)^2
       =\sum_{m\ge1}a_mq^m.
\]
Shell counts, Fourier coefficients and Hecke operators at nonintegral
positive indices are interpreted as zero. The full signed ordered shell
count satisfies
\[
 r_F(m)=\frac{12}{5}\{\sigma(m)-11\sigma(m/11)\}
                       +\frac85a_m\quad(m\ge1).           \tag{M11.4}
\]
Here \(f_{11}\) is the normalized weight-two newform of level 11
(Martin--Ono, Theorem 1, level-11 row), so for \(11\nmid m\),
\[
                         |a_m|\le\tau(m)\sqrt m.           \tag{M11.5}
\]
This is the unramified weight-two Ramanujan bound and the all-ideal
prime-power recurrence, with the same exact imports as the spectral
argument below.

Here is the complete finite proof of (M11.4). The theta series of the
rank-four even Gram lattice has level 11 and trivial character. The
Eisenstein series
\((11E_2(11\tau)-E_2(\tau))/10\) is holomorphic of weight two at that
level. Indeed it equals
\[
 \frac1{10}\frac1{2\pi i}\frac{d}{d\tau}
                     \log\frac{\Delta(11\tau)}{\Delta(\tau)},
 \qquad \Delta=\eta^{24}.
\]
The quotient is a modular function, nonvanishing in the upper
half-plane and meromorphic at the cusps. In a cusp coordinate its
expansion is \(q_c^\nu U(q_c)\), \(U(0)\ne0\); its logarithmic
derivative is a constant plus \(q_c U'/U\), up to the fixed cusp-width
factor, so it has no cusp pole. Subtract this Eisenstein series and
\(8f_{11}/5\) from the theta series and multiply by 5. The result has
integer coefficients: the Eisenstein series multiplied by 5 has
constant 5 and positive coefficients
\(12\{\sigma(n)-11\sigma(n/11)\}\). The Sturm bound is
\((2/12)[\mathrm{SL}_2(\mathbb Z):\Gamma_0(11)]=2\).
The theta coefficients through \(q^2\) are \((1,4,4)\), and the cusp
coefficients are \((0,1,-2)\). The constant and the first two
coefficients of the difference vanish. Sturm's congruence theorem
applied modulo every rational prime proves equality. The checker uses
independent lattice enumeration and truncated eta-product multiplication
for this complete certificate, and challenges it further through 256.

Because \(R/2R=\mathbb F_4\), an even field norm has both coordinates
even. For an even target \(n\), either both norms are odd or all four
coordinates are even. Thus the count \(G(n)\) of shell points with both
norms odd is exactly
\[
 G(n)=r_F(n)-r_F(n/4).
\]
Put
\[
 S(n)=\sigma(n)-\sigma(n/4).
\]
For an even \(n=2^t m\) with \(m\) odd and \(11\nmid n\), define
\[
 \Delta(n)=\tau(m)\sqrt m\,2^{t/2}\frac{t+3}{2},\qquad
 E=\frac{2\Delta(n)}{3S(n)}.
\]
For every such target,
\[
 S(n)=3\cdot2^{t-1}\sigma(m)\ge\frac32 n,\qquad
 G(n)\ge\frac{12}{5}S(n)(1-E).                            \tag{M11.6}
\]
The bound on the filtered cusp coefficient used here is
\(|a_n-a_{n/4}|\le\Delta(n)\). To see this, let the two normalized
Hecke roots at 2 be \(u,v\), with \(uv=1\) and \(|u|=|v|=1\).
The all-ideal recurrence gives
\[
 a_{2^t}-a_{2^{t-2}}
 =2^{t/2}\left(u^t+v^t+
       \frac12\sum_{j=1}^{t-1}u^jv^{t-j}\right).
\]
The second coefficient is zero for \(t=1\), where the displayed
sum is empty. For \(t\ge2\), the subtraction cancels half of every
interior term because \(uv=1\). The triangle inequality gives the
factor \((t+3)/2\); multiplicativity supplies \(\tau(m)\sqrt m\).
This proves (M11.6) and handles arbitrarily deep two-adic valuations.
The restriction \(11\nmid n\) is necessary: at \(n=22\), direct
enumeration gives \(G(22)=4\), while the displayed lower bound would be
\(432/5-(32/5)\sqrt{22}>56\). The single permitted ramified factor is
removed before applying this estimate and restored only to the final
value pair.
The coarser bound \(E\le\tau(n)\sqrt n/S(n)\) also follows from
\(\Delta(n)\le(3/2)\tau(n)\sqrt n\).

## Hecke operators with all ideal classes retained

The finite quotient must retain every locally principal right-ideal
class of \(O\). The shell itself is a return to the principal class,
rather than the entire Hecke row. We give the normalization explicitly.

For projective unit stabilizers \(\Gamma_i\) of the maximal left orders,
Voight's Theorem 25.1.1 gives
\[
 \mathfrak m=\sum_i\frac1{|\Gamma_i|}=\frac{11-1}{12}=\frac56.
\]
The principal class has \(|\Gamma_0|=2\), hence probability
\[
                         w=\frac{1/2}{5/6}=\frac35.        \tag{M11.7}
\]
No assertion that \(O\) has class number one is used.

At an odd unramified prime \(p\), reduce projective maximal-order units
to \(G_p=\mathrm{PGL}_2(\mathbb F_p)\). At 11 use the radical quotient
\[
 G_{11}=\mathbb F_{121}^{\times}/\mathbb F_{11}^{\times},
                           \qquad |G_{11}|=12.
\]
Let \(X_p\) be the adelic quotient with the reduction kernel imposed
at \(p\), and maximal unit level everywhere else. Its fibers over
ideal classes are \(\Gamma_i\backslash G_p\); at 11 use \(G_{11}\).
Every \(\Gamma_i\) injects into the specified finite group. Indeed a
noncentral norm-one unit has integral trace \(-1,0,1\), so
\(\operatorname{tr}(u)^2-4=-3\) or \(-4\); scalar reduction would
force divisibility by \(p\) (by \(p^2\) at a split prime). Our odd
inert primes are at least 7, and the ramified prime is 11.
All fine stabilizers are consequently trivial. Haar probability is
uniform counting probability on \(X_p\), with
\[
 |X_p|=|G_p|\mathfrak m,\qquad
 |\Gamma_0\backslash G_p|/|X_p|=w.
\]
In particular \(X_{11}\) has ten points and its principal fiber has six.

Let \(T_m\) be the all-integral-ideal Hecke operator for
\((m,11p)=1\), or \((m,11)=1\) on \(X_{11}\).
It is self-adjoint for this invariant probability and commutes with the
right finite group action. If \(h\) is supported in the principal fiber,
then at its base point \(x_0\),
\[
 (T_mh)(x_0)=\frac14\sum_{\substack{\alpha\in O\\\operatorname{nrd}\alpha=m}}
                              h(\rho_p(\alpha)).           \tag{M11.8}
\]
Only principal-return ideals contribute. Generators differ by four
units; principalization initially gives \(\alpha^{-1}\), but
projectively \(\alpha^{-1}=\bar\alpha\) and conjugation permutes the
shell. Hence the stated orientation and factor agree.

The constant eigenvalue is \(\sigma(m)\); the quadratic norm character
\(\psi_p\) has eigenvalue \((m/p)\sigma(m)\). At 11 the latter symbol
means \((m/11)\). These are the only global norm characters. One way to
check this without a class-one hypothesis is to pass to reduced norms:
away from \(p\) the maximal units give \(\mathbb Z_\ell^\times\), at
\(p\) the reduction kernel gives square residues, central scalars give
squares, and rational norms give all positive rational numbers by the
Hasse--Schilling theorem (Voight, Main Theorem 14.7.4). The remaining
norm quotient has order two. Strong approximation for the
norm-one group at a split place excludes extra constituents that are
characters at a split place.

On the orthogonal complement of \(1,\psi_p\),
\[
                         \|T_m\|\le\tau(m)\sqrt m.         \tag{M11.9}
\]
To specify the imported argument, the proof of the
[Eisenstein spectral lemma](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex)
uses the all-ideal correspondence, strong approximation, global
Jacquet--Langlands, and the weight-two Ramanujan bound. The same proof
applies on this multi-class quotient: retaining all fibers replaces
class-one transitivity by Haar measure and (M11.7); the scalar kernel
was just checked; the only character constituents were just identified;
the trivial infinite type transfers to weight two; every prime dividing
\(m\) is split in \(B\) and has maximal local level. At each such prime
\(\ell\), the all-ideal eigenvalue at \(\ell^a\) is the sum of \(a+1\)
monomials in Hecke roots of absolute value \(\sqrt\ell\). Multiplication
over primes gives (M11.9). The exact imported theorem editions and the
existing proof bytes are pinned in the source ledger. The bound is not
being asserted at an added level prime dividing \(m\).

On any of these quotients the odd-component filter is
\(D_n=T_n-T_{n/4}\), with the second operator zero for \(4\nmid n\).
Scalar multiplication by 2 is projectively trivial, so (M11.8) identifies
this difference with the filtered counts. Its constant and character
eigenvalues are \(S(n)\) and \((n/p)S(n)\), and its complementary norm is at most \(\Delta(n)\). The same
Satake-root cancellation just proved applies to each simultaneous cusp
eigenvector, and self-adjointness gives the operator norm.
There is also a direct spectral derivation of the denominator in
(M11.6). On the coarse ideal-class space the principal indicator has
mean \(w=3/5\). Its nonconstant squared norm is \(w(1-w)\), and
only the constant norm character occurs at full maximal-unit level.
Evaluating its \(D_n\) correlation at the principal class gives
\(G(n)/4\ge wS(n)-(1-w)\Delta(n)\), exactly (M11.6).
This derivation independently checks the theta-based normalization.


## Small bad-prime counts and the ramified filter

At an odd inert \(p\nmid n\), the normalizer \(H_p\) of the nonsplit
torus from \(R/pR\) has order \(2(p+1)\), contains \(\Gamma_0\), and
in the principal fiber describes exactly those invertible quaternions
with a zero field component. Put
\[
 \mu_p=\frac{2}{p(p-1)},\quad
 h=1_{\Gamma_0\backslash H_p},\quad \nu=w\mu_p.
\]
The mean of \(h\) is \(\nu\). The determinant character has mean zero
on this orbit, since the torus norm is onto \(\mathbb F_p^\times\).
Both \(h\) and \(D_nh\) are right-\(H_p\)-invariant. Averaging at the
base point therefore gives
\[
 (D_nh)(x_0)=\frac{\langle h,D_nh\rangle}{\nu}.
\]
The complementary squared norm is \(\nu(1-\nu)\). By (M11.8)--(M11.9),
the bad filtered count satisfies
\[
 \frac{b_p^{\rm odd}(n)}4
       \le w\mu_p S(n)+\Delta(n),
 \qquad
 \frac{b_p^{\rm odd}(n)}{G(n)}
       \le\frac{\mu_p+(5/2)E}{1-E}.                        \tag{M11.10}
\]
If \(p\mid n\), anisotropy forces both field components divisible by
\(p\) as soon as one is. The count is zero unless \(p^2\mid n\), and
otherwise equals \(G(n/p^2)\). From (M11.4), its ratio is at most
\((p^{-2}+E/p)/(1-E)\), since
\(S(n/p^2)/S(n)\le p^{-2}\) and
\(\Delta(n/p^2)\le\Delta(n)/p\). This is bounded by (M11.10), covering all
target valuations without an operator at its own level prime.

At 11 the residue order is \(\mathbb F_{11}[j]=\mathbb F_{121}\).
The two field-component axes project to \(\{1,j\}=\Gamma_0\) in
\(G_{11}\). The bad subset of the principal fiber is thus the single
point \(x_0\). Its indicator \(\delta\) has mean \(1/10\), and its
projections on \(1\) and \(\psi_{11}\) both have coefficient \(1/10\).
Its complementary squared norm is \(1/10-2/100=2/25\).
Evaluating the correlation at a point of probability \(1/10\) gives
\[
 \frac{b_{11}^{\rm odd}(n)}4
    \le\frac{1+(n/11)}{10}S(n)+\frac45\Delta(n),
 \qquad
 \frac{b_{11}^{\rm odd}(n)}{G(n)}
    \le\frac{1/3+2E}{1-E}.                                \tag{M11.11}
\]
The worst character sign was used in the ratio.

## Large primes and the positive sieve margin

Define
\[
 M(X)=\max_{1\le t\le X}
              \prod_{\substack{q\mid t\\(q/11)=1}}(v_q(t)+1),
 \qquad \mathcal I=\{p\text{ odd prime}:(p/11)=-1\}.
\]
Class number one and ideal factorization give \(r_Q(t)\le2M(X)\) for
\(1\le t\le X\). Also there are at most \(6X\) nonzero \(R\)-vectors
of norm at most \(X\). Here is an elementary bound with that constant.
The form dominates \(Q_E=x^2+xy+y^2\). Each six-unit orbit for \(Q_E\)
has exactly one representative \(x\ge1,y\ge0\). For \(s=x+y\), assign
it the distinct integer rank \(s(s-1)/2+x\). For \(s\ge2\),
\[
 s(s-1)/2+x\le s(s+1)/2\le\lceil3s^2/4\rceil\le Q_E(x,y);
\]
\(s=1\) is immediate. At most \(\lfloor X\rfloor\) ranks occur.

Both field components in the odd shell are nonzero. If an inert prime
\(p\) divides a component norm, that component is \(p\) times a nonzero
vector, giving at most \(6n/p^2\) choices. There are at most \(2M(n)\)
choices for the other component. Counting both orientations gives
\[
 b_p^{\rm odd}(n)\le24M(n)n/p^2.                           \tag{M11.12}
\]
For \(P>0\), write
\[
 B=\sum_{p\in\mathcal I}\mu_p,\quad
 \Pi(P)=\#\{p\in\mathcal I:p\le P\},\quad
 T(P)=\sum_{\substack{p\in\mathcal I\\p>P}}p^{-2}.
\]
The union bound, using (M11.6), (M11.10)--(M11.12), leaves at least
\[
 \frac{12}{5}S(n)
 \left\{\frac23-B-\frac{20}{3}M(n)T(P)
                         -E\left(3+\frac52\Pi(P)\right)\right\}
                                                               \tag{M11.13}
\]
points. In deriving this, the denominators \(1-E\) cancel. The final
bounds below ensure \(E<1\). Positive surviving norms belong to \(C\).
Split common factors in these coordinate witnesses cause no problem:
(M11.3) supplies primitive representations of their values separately.

## Exact constants and the two cutoffs

Odd inert primes occupy the five residue classes
\(7,13,17,19,21\pmod{22}\). Counting these progressions and integrating
\(x^{-2}\) on each gives
\[
 \Pi(P)\le\frac5{22}P+\frac32,\qquad
 T(P)\le\frac5{22P}+\frac5{P^2},\qquad B<\frac9{100}.       \tag{M11.14}
\]
For the last inequality, sum the actual inert primes through 64 and
majorize the rest by the five progressions starting at
\(65,73,79,83,85\). Since \(p\ge7\),
\(\mu_p\le(7/3)p^{-2}\); each progression beginning at \(a\) contributes
at most \(a^{-2}+(22a)^{-1}\). All checks are exact rational arithmetic.

The finite Euler certificates prove, for all real \(X\ge1\),
\[
 M(X)\le7X^{1/6},\qquad
 \tau(m)\le29m^{1/6}\ (m\text{ odd},11\nmid m),\qquad
 \tau(n)\le15n^{1/6}\ (p\mid n\Rightarrow(p/11)=-1).        \tag{M11.15}
\]
For each relevant prime \(p<64\), maximize \((a+1)^6/p^a\) for
\(0\le a\le10\), then multiply these rational maxima. Larger exponents
cannot improve them, because \(12^6<2\cdot11^6\) and the successive
ratio decreases with \(a\). For \(p\ge64\), \(a+1\le2^a\) makes the
factor at most one. The three products are respectively below
\(7^6,29^6,15^6\). This is the completeness argument for the finite
certificates and the unrestricted bounds (M11.15).

For arbitrary even \(n=2^t m\) prime to 11, (M11.6) and the odd
bound in (M11.15) give
\[
 E\le\frac{58}{9}(t+3)2^{-t/6}n^{-1/3}\le29n^{-1/3}.
\]
The last maximum occurs at \(t=6\), where
\((t+3)2^{-t/6}=9/2\). This holds for every \(t\ge1\): the consecutive
ratios decrease with \(t\), and
\(9^6>2\cdot8^6\), \(10^6<2\cdot9^6\) bracket the change from
increase to decrease. Choose \(P=(4/5)n^{1/4}\). The terms after
\(B\) in (M11.13) are at most
\[
 \frac{1745}{66}n^{-1/12}+\frac{1681}{3}n^{-1/3}.
\]
For \(n\ge2^{67}\), use \(n^{1/12}>47\), since \(47^{12}<2^{67}\),
and \(n^{1/3}>2^{22}\). The resulting exact bound is
\[
 \frac9{100}+\frac{1745}{66\cdot47}+\frac{1681}{3\cdot2^{22}}
 =\frac{106147055309}{162634137600}<\frac23.
\]
If \(N\ge2^{71}\), then \(n=N/11^e>2^{67}\). Thus (M11.13) is
positive and proves the existence part of Theorem M11.
Its margin exceeds \(3/1000\). The number of surviving coordinate
points is at least \((18/5)n(3/1000)\), while a fixed value pair has
at most \((2M(n))^2\le196n^{1/3}\) witnesses. There are therefore at
least \(n^{2/3}/20000\) ordered value pairs. Multiplication by \(11^e\)
is injective on these pairs, and \(11^{2/3}<5\), proving (M11.1).

If \(n\) has only inert prime factors, the last bound in (M11.15) gives
\(E\le10n^{-1/3}\). With \(P=n^{1/4}\), the criterion terms after \(B\)
are bounded by
\[
 \frac{1075}{66}n^{-1/12}+\frac{1805}{6}n^{-1/3}.
\]
For \(n\ge2^{60}\),
\[
 \frac9{100}+\frac{1075}{66\cdot32}+\frac{1805}{6\cdot2^{20}}
       =\frac{345616637}{576716800}<\frac23.               \tag{M11.16}
\]
Finally multiplication of a member of \(P\) by any split prime preserves
membership, by (M11.3). If a least admissible counterexample had a split
prime factor \(q\), then \(N/q\) would still be even, at least 2, and
not divisible by 121. A representation of \(N/q\) would lift to one of
\(N\), a contradiction. A least counterexample has only inert factors
and at most one factor 11; (M11.16) forces its reduced part below
\(2^{60}\). Conversely every failure in (M11.2) is a failure of the
identity. This proves the intermediate finite equivalence. The following section
exhausts it by analytic intervals and a complete bounded certificate.

## Closing the finite core: sharp bounds and a complete certificate

The sharpness needed here comes from the binary lattice count and the
actual factor \(\sigma(m)\) in the denominator. No stronger automorphic
input is required.

For every real \(X\ge0\),
\[
 \#\{(x,y)\in\mathbb Z^2:0<Q(x,y)\le X\}\le\frac{12}{5}X.
                                                               \tag{M11.17}
\]
To prove this, put \(g(y)=\sqrt{\max(X-11y^2/4,0)}\). Each integer
row contains at most \(2g(y)+1\) points. Its integral and monotonicity
give, after deleting the origin,
\[
 S_Q(X)\le\frac{2\pi}{\sqrt{11}}X+
                         (2+4/\sqrt{11})\sqrt X.
\]
For integer \(X\ge49\), use \(\pi<22/7\), \(\sqrt{11}>33/10\), and
\(\sqrt X\ge7\): the resulting ratio is less than
\(440/231+106/231=26/11<12/5\). Exact lattice enumeration verifies
\(5S_Q(j)\le12j\) for each \(1\le j\le48\); equality holds at 5.
For general real \(X\), apply the integer bound at \(\lfloor X\rfloor\).
Thus the coefficient \(20/3\) in (M11.13) improves to \(8/3\): the
large-prime numerator is now at most \((48/5)M(n)n/p^2\).

For an even inert core \(n=2^t m\), \(m\) odd, the exact definition
of \(E\) gives
\[
 E n^{1/3}=
 \frac{2(t+3)}9\,2^{-t/6}
 \prod_{p^a\parallel m}
                  \frac{(a+1)p^{5a/6}}{\sigma(p^a)}
 \ <\ \frac{11}{4}.                                      \tag{M11.18}
\]
The dyadic factor is at most one, with its maximum at \(t=6\), as
proved above. For odd inert \(p\ge64\), the sixth power of the Euler
factor is at most \((a+1)^6/p^a\le1\). For odd inert \(p<64\),
the ratio of consecutive sixth powers, once \(a\ge2\), is strictly
less than \((4/3)^6/p\le(4/3)^6/7<1\), because
\(\sigma(p^{a+1})>p\sigma(p^a)\). Hence only exponents 0,1,2 at
\(p=7,13,17,19,29,41,43,61\) need be maximized. The exact product is
\[
 \frac{184167468895316154074204818561020107}
      {465817612584447478971000000000000}
                 <(11/4)^6.
\]
This proves (M11.18) for every inert core, without a bound on its size.
It is not asserted for targets containing split prime factors.

The envelope \(M(X)\) is also exactly computable on the bounded
intervals needed here. Remove nonsplit prime factors from a maximizing
integer. Assign its nonzero exponents in decreasing order to the
smallest split primes. This preserves its divisor product and decreases
its size: exchanging exponents \(a<b\) at primes \(p<q\) multiplies
the integer by \((p/q)^{b-a}<1\). Therefore it suffices to enumerate
initial split-prime prefixes with nonincreasing positive exponents.
The recursive enumeration checks the size bound before each extension.
The finite prime and exponent lists are complete because the first 12
split primes already have product greater than \(2^{61}\), and
\(3^{41}>2^{61}\). The retained checker verifies these inequalities.
Sorting the resulting size/divisor-product pairs gives an exact maximum
for every argument below \(2^{60}\).

For each \(k=25,\ldots,59\), put
\[
 L_k=2^k,\quad U_k=2^{k+1}-1,\quad d_k=M(U_k),\quad
 e_k=\frac{11}{4\lfloor\sqrt[3]{L_k}\rfloor},\qquad
 P_k=\left\lfloor\sqrt{\frac{16d_k}{15e_k}}\right\rfloor.
\]
Thus throughout that interval \(M(n)\le d_k\) and \(E\le e_k<1\).
Let \(R=65536\), \(D=2^{64}\). Sieve all primes through \(R\)
exactly, retaining the odd inert primes, and set
\[
 A_k=\#\{p\in\mathcal I:p\le P_k\},\quad
 V_k=\frac1D\sum_{\substack{p\in\mathcal I\\P_k<p\le R}}
                   \left\lceil\frac{D}{p^2}\right\rceil
                             +\frac5{22R}+\frac5{R^2}.
\]
Every \(P_k<R\), so \(A_k=\Pi(P_k)\) and \(V_k\ge T(P_k)\).
The final two terms bound every prime beyond the finite sieve using
(M11.14). The finite certificate verifies, in exact rational arithmetic,
all 35 inequalities
\[
 \frac9{100}+\frac83d_kV_k+e_k\left(3+\frac52A_k\right)
                       <\frac23-\frac1{1024}.             \tag{M11.19}
\]
The worst interval is \(k=25\); all individual bounds and inputs are
retained in the receipt. Consequently (M11.13), with the improved
coefficient, proves \(n\in C+C\) for every inert core
\(2^{25}\le n<2^{60}\). Formula (M11.16) covers every inert core at
least \(2^{60}\); it also gives \(C+C\), since the sifted norms are
units at 11. These intervals leave no unbounded or untested gap.

It remains to cover the small interval. The checker constructs the
exact flags of \(C\cap[1,2^{25}]\) by starting with all positive
integers and deleting multiples of every forbidden prime, including 2
and 11. This is the semigroup definition, not a heuristic primality
condition. Independent primitive coordinate enumeration, with 11-units
retained, agrees with these flags through 10000.

Pack the flags into a binary integer \(B_C\), with bit \(j\) indicating
\(j\in C\). The exact union
\[
 \bigvee_{\substack{a\le5000\\a\in C}} (B_C\ll a)          \tag{M11.20}
\]
uses 774 short summands. After truncation through \(2^{25}\), every
one of the
\[
 2^{24}-\lfloor2^{25}/22\rfloor=15252015
\]
eligible even, 11-unit target bits is present. Thus every eligible
target in this full interval has two positive summands in \(C\).
No claim that a short summand must exist was assumed: the computed
union covers every required bit. The receipt retains the exact domain,
summand count, zero-missing outcome and flag/coverage hashes; the checker
reconstructs the flags and union. The hash alone is not the proof.

Finally suppose a least even 11-unit exception to \(C+C\) existed.
It exceeds \(2^{25}\) by (M11.20). If it had a split prime factor
\(q\), dividing by this odd prime would leave a smaller eligible even
target at least 2. Multiplying its two \(C\)-summands by \(q\) would
lift a representation, since \(C\) is multiplicatively closed. The least
exception must therefore be an inert core, which (M11.19) and (M11.16)
exclude. Conversely \(C+C\) has no odd target, and no multiple of 11:
its two residues are nonzero squares and \(-1\) is nonsquare modulo 11.
This proves the first identity (M11.0).

For an even \(N\) with \(121\nmid N\), write \(N=11^e n\), \(e=0,1\).
Then \(n\in C+C\), and multiplying both summands by \(11^e\) puts
them in \(P=C\cup11C\). Odd targets and multiples of 121 were already
excluded locally. This proves the second identity and completes the
exact classification. \(\square\)

## Evidence and scope

The [checker](verify_primitive_norm_minus11.py) and its
[receipt](primitive_norm_minus11_checks.json) certify the complete Sturm
coefficient list, finite Euler maxima, sharp lattice bounds, all 35
rational interval inequalities and the exhaustive small-core coverage.
They also challenge the theta identity and ramified coordinate map by
bounded direct enumeration. The [source ledger](primitive_norm_minus11_source_ledger.json)
pins Voight's actual book edition, the primary eta/Sturm sources and the
existing spectral proof whose imported theorem statements are reused.
Run from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_minus11.py \
  --output research/arithmetic/primitive_norm_minus11_checks.json
```

Certificate verification requires active Python assertions. The checker rejects optimized execution (-O, -OO, or optimization enabled through PYTHONOPTIMIZE) before generating a receipt.

The separate [coordinate membership verifier](verify_primitive_norm_coordinate_membership.py)
authenticates every core flag through \(2^{25}\). Simultaneous sign change
allows \(y\ge0\), and the reflection \(x\mapsto-x-y\) allows
\(2x+y\ge0\), preserving the norm and gcd. The identity
\(4Q=(2x+y)^2+11y^2\) gives exact enclosing integer limits.
Keeping positive odd 11-unit norms with gcd one therefore enumerates
the complete core independently of prime-support sieving. The two
full bytearrays are compared directly before hashing; coverage-mask
equality alone would not authenticate the summands. The shared
[membership receipt](primitive_norm_coordinate_membership_checks.json)
records this check and rejection of false high membership flags.
Run from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_coordinate_membership.py \
  --output research/arithmetic/primitive_norm_coordinate_membership_checks.json
```

The infinite spectral estimate is a written deduction from those imports,
not a consequence of the finite checker passing. Internal agents checked
the local and global normalizations, exact interval coverage and descent
separately, and replayed the small-core computation. This is internal
review, without formal verification, external review or an established
priority claim. The separate [exact -8 proof](PRIMITIVE_NORM_MINUS8_EXACT.md) now settles
that binary identity too. The independent \(C_{12}+C_{12}\) problem remains open.
