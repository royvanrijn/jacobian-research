# Exact core sumsets at five Heegner discriminants

Let \(p\in\{7,19,43,67,163\}\). Let \(C_p\) consist of 1 and positive
integers whose prime divisors \(q\) satisfy \((q/p)=1\).
Let \(L_p\) consist of the integers \(n\ge2\) with \(p\nmid n\),
requiring \(n\) even exactly when 2 is inert in the field.

**Theorem.** For all five primes, \(C_p+C_p=L_p\setminus E_p\),
with the complete exception sets below. There is a common effective
sieve bound \(2^{120}\), independent of the finite scans. The exact
lists give the sharper common threshold 993 for eligible targets.

| \(p\) | Core exceptions | Largest | Exceptions with a split divisor |
| --- | ---: | ---: | --- |
| 7 | \(\{11,29\}\) | 29 | 11, 29 |
| 19 | \(\{4,20\}\) | 20 | \(4\cdot5\) |
| 43 | 10 | 74 | none |
| 67 | 21 | 138 | \(4\cdot17,\ 6\cdot19,\ 4\cdot29,\ 6\cdot23\) |
| 163 | 67 | 992 | \(4\cdot41,\ 4\cdot47,\ 4\cdot83,\ 8\cdot43,\ 8\cdot71\) |

The full lists at \(p>7\) are in the [exception receipt](heegner_core_sumset_checks.json).
Here \(E_{43}=\{4,6,8,10,16,20,38,50,56,74\}\) and
\[
 E_{67}=\{4,6,8,10,12,14,16,22,26,28,32,44,50,62,68,80,86,98,114,116,138\}.
\]
At \(p=163\), the original scan through \(2^{34}\) is joined to the
proved infinite bounds by the factorization search below.
[Completion receipt](heegner_core_completion_checks.json)

The [prime-discriminant extension](PRIME_DISCRIMINANT_CORE_UNIFORM.md)
proves a uniform polynomial bound and the limiting finite-window gap law.

## Comparison with discriminant −7

Write \(P_p\) for the positive primitive values of
\(Q_p(x,y)=x^2+xy+(p+1)y^2/4\). The quadratic orders here have class number
one. For example, the imaginary quadratic Minkowski bound
\(2\sqrt p/\pi\) leaves only inert rational primes below the bound in each
of the four cases. Thus ideal factorization gives
\(P_p=C_p\sqcup pC_p\). At a split prime an oriented principal factor can
be raised to any power. At the ramified prime one can use
\(\sqrt{-p}\) once. An inert prime dividing a norm divides both
coordinates.

For \(p>7\) here, 2 is inert, so core and primitive values are odd.
At \(p=7\), 2 splits and there is no dyadic restriction. The core
values are nonzero squares modulo \(p\), and \(-1\) is a nonsquare.
Conversely, for \(n\ne0\bmod p\) there are
\(p+1-2(1+(n/p))>0\) ordered solutions to \(x^2+y^2=n\bmod p\) with both
variables nonzero. At any other odd prime local unit conditions on both summands exclude
at most two residue choices. These choices lift. For \(p=7\) and an
odd target, take one summand twice an odd value and the other odd.
The factor 2 is split and multiplication by 2 preserves the square
classes modulo 7. After dividing that even summand's prescribed
residues by 2, the Chinese remainder theorem and Dirichlet's theorem
realize the joint unit residues by split primes. For even targets
both summands can be odd. Thus \(L_p\) is the full congruence-admissible
set.

If \(p\mid N\), the same nonsquare argument forces both primitive
summands to be divisible by \(p\), giving the exact identity
\[
 (P_p+P_p)\cap p\mathbb Z=p(C_p+C_p).
\]
For \(p\nmid N\), the two possibilities are \(C_p+C_p\) and \(C_p+pC_p\).
Consequently the exceptions among even \(N\ge2\) with \(p^2\nmid N\) are
\[
 p E_p \ \sqcup\ \bigl(E_p\setminus(C_p+pC_p)\bigr)
\]
when \(E_p\) is the full core exception set. The corresponding bounded
identity uses the core scan up to the chosen target bound and dilation
only up to that bound.

At \(p=7\) the same dilation identity gives \(7\{1,11,29\}\), with 1
included because the primitive target 7 is locally admissible and the
summands are positive. Both core exceptions 11 and 29 are split primes,
and the mixed sums \(4+7\) and \(1+28\) fill their undilated gaps.
The four cases \(p>7\) have 1, 8, 17 and 62 undilated exceptions.
For instance, \(20=1+19\) fills one core gap at \(p=19\), while 4 survives.
The complete primitive exception counts at \(p>7\) are 3, 18, 38 and 129. At \(p=19\) the full primitive list is
\(\{4,76,380\}\).
The \(p=11\) comparison has the same core congruence condition as the
four new cases and no core gaps.
[Discriminant −7](PRIMITIVE_NORM_MINUS7_EXACT.md) ·
[Discriminant −11](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md)

The least core value greater than 1 is \((p+1)/4\).
Indeed \(4Q_p(x,y)=(2x+y)^2+py^2\), the primitive \(y=0\) value is 1,
and \((x,y)=(0,1)\) attains the bound. This forces every even
\(4\le n<1+(p+1)/4\) to be a gap. Beyond that initial range the lists
depend on which values are split-supported.
If \(q\) splits and \(q\mid n\in E_p\), then \(n/q\in E_p\) because
multiplying a representation by \(q\) lifts it. This divisibility
property holds throughout the exact scan and also for the full set.
In the observed lists, every exception with a split divisor has
exactly one such prime factor. Its quotient is always 4, 6 or 8.


## Split divisibility and the baseline for inert gaps

Write an admissible target \(n\) as \(r c\), with \(r\) supported on
inert primes and \(c\in C_p\). For \(p>7\), \(r\) is even.
For each fixed \(r\), the set of \(c\) for which \(rc\) is a gap
is closed under taking divisors in \(C_p\). Indeed any representation
at \(rc'\) with \(c'\mid c\) scales to one at \(rc\).
In particular a split-bearing gap must lie above an inert-supported
gap. At \(p=7\) this holds after adjoining 1 to the gaps, which
accounts for the positivity boundary and permits \(r=1\).

There is a second constraint on a prime multiple \(rq\) of a gap \(r\).
If \(rq=a+b\) is represented and \(q\mid a\), then \(q\mid b\) and
division by \(q\) would represent \(r\). Thus any representation of
\(rq\) must have both summands prime to \(q\). The ten failures of this
type at \(p>7\) occur in the following fibers.

| \(p\) | Inert part \(r\) | Split primes \(q\) with \(rq\) a gap |
| --- | ---: | --- |
| 19 | 4 | 5 |
| 67 | 4 | 17, 29 |
| 67 | 6 | 19, 23 |
| 163 | 4 | 41, 47, 83 |
| 163 | 8 | 43, 71 |

For each failure, the completion receipt lists every possible shorter
core summand and a forbidden prime divisor of its complement.
For example, with \(p=19\) the only possibilities below or at \(20/2\)
are \(1,5,7\). Their complements \(19,15,13\) are blocked by the
ramified prime 19 or the inert primes 3 and 13. This explains the
individual gaps using additive information beyond split divisibility.

To quantify the observed 90 inert-supported gaps among the 100 gaps
for \(p>7\), take for each \(p\) the window up to its largest exception,
restrict to admissible targets, and condition on the observed number
of gaps for that \(p\). Uniform sampling of targets within those
windows gives an expected 77.52 inert-supported gaps. That model
allows split-bearing gaps without their inert roots.

A second model samples uniformly among divisor-closed sets in the
same windows and of the same sizes, requiring the observed last gap
and all decisions forced below twice the least nonunit core value.
In that initial range \(n=1+(n-1)\) is the only possible positive
core decomposition up to order. Every fiber in the four windows
has at most one split prime per target, counted with multiplicity.
Thus it consists of an inert root and its available prime multiples.
A chosen root with \(f\) required children and \(d\) optional children
contributes \(y x^{1+f}(1+x)^d\) to the counting polynomial.
The variable \(x\) counts gaps and \(y\) counts inert roots.
The empty-fiber term is allowed when neither the root nor any child
is required. A forced represented root permits only the empty state.
Extract the coefficient at the observed number of gaps in each
window, and combine the four distributions. This model has expected
inert count 90.27 and probability about \(0.6383\) of at least 90.
These probabilities describe the specified model of divisor-closed
sets. They do not assign probabilities to the arithmetic gap sets.
[Exact baseline and local nonrepresentation certificates](heegner_core_completion_checks.json)

## Quaternionic sieve

The argument uses the all-ideal Hecke normalization in the −11 proof.
Here \(O=\mathbb Z[(1+\sqrt{-p})/2]\oplus
\mathbb Z[(1+\sqrt{-p})/2]j\) in \((-p,-1)_{\mathbb Q}\), where
\(j^2=-1\) and \(jz=\bar z j\). It is a maximal order of discriminant
\(p\). Its norm is \(Q_p\oplus Q_p\) and its unit group is
\(\{1,-1,j,-j\}\), by the minimum of \(Q_p\) just established.
The mass formula over all right-ideal classes gives principal-class
probability \(w=6/(p-1)\).
[Voight, Example 15.5.7 and Theorem 25.1.1](sources/voight-quaternion-algebras-author.pdf)

At an odd auxiliary prime \(q\ne p\) use the projective maximal-order
unit quotient \({\rm PGL}_2(\mathbb F_q)\). At \(p\) use
\(\mathbb F_{p^2}^{\times}/\mathbb F_p^{\times}\) from the radical
quotient. The fiber over each ideal class is the corresponding
projective unit group quotient. Its stabilizer injects into the finite
group. A noncentral unit has trace discriminant \(-3\) or \(-4\).
Scalar reduction at an odd split place of the quaternion algebra
would force divisibility of the discriminant by \(q^2\), which excludes
even \(q=3\). At \(p\) scalar reduction would force divisibility by
\(p\), also impossible. Thus the fine quotients carry uniform counting
probability and the coarse quotient carries inverse-unit mass.

As in the −11 proof, the all-ideal Hecke operator \(T_m\) is
self-adjoint and has principal-return shell sum divided by four.
It commutes with the right level action. Away from the discriminant
and the auxiliary level, its eigenvalue on constants is \(\sigma(m)\).
At one auxiliary level the only extra norm character is the quadratic
character at that level. Its eigenvalue is \((m/q)\sigma(m)\), using
\(q=p\) for the radical level. These are the norm characters because
the level condition restricts the local unit norms to square residue
classes at that one prime, and Hasse–Schilling and strong approximation
give the norm quotient as in the −11 proof. On the complement the
Jacquet–Langlands transfer at trivial infinite type and the weight-two
Ramanujan bound give operator norm at most \(\tau(m)\sqrt m\).
The transfer uses maximal local level at every prime dividing the
Hecke index, including 2. The estimates below use separate shell
reduction when the auxiliary prime divides the target.
[Multi-class operator normalization](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md#hecke-operators-with-all-ideal-classes-retained) ·
[Classical Hecke inputs](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex)

For \(n=2^t m\in L_p\), with \(m\) odd, apply the odd-component filter
\(D_n=T_n-T_{n/4}\), with the second term zero at a noninteger index.
At 2 the field norm is anisotropic, so in an even shell both blocks
are odd or both are divisible by 2 as vectors. The constant eigenvalue
and a bound on complementary eigenvalues are
\[
 S=3\cdot2^{t-1}\sigma(m),\qquad
 \Delta=2^{t/2}(t+3)\tau(m)\sqrt m/2 .
\]
The \((t+3)/2\) comes from subtracting the two all-ideal Satake sums
at 2, as in the −11 proof. Set \(E=\Delta/S\). The principal-class
indicator gives an odd shell with at least
\(4[wS-(1-w)\Delta]\) signed ordered points.

For an odd inert \(q\nmid n\), a bad point has one zero block modulo
\(q\). Such invertible elements project to the normalizer of the
nonsplit torus. Its indicator on the principal fiber has mean
\(w\mu_q\), where \(\mu_q=2/(q(q-1))\), and zero quadratic norm character
component. Average the evaluation at the base point over that
normalizer. The resulting self-correlation bounds the bad count
divided by four by \(w\mu_q S+(1-w\mu_q)\Delta\).
If \(q\mid n\), a bad point forces both blocks divisible by \(q\).
The count is zero at valuation one. At valuation at least two it
is bounded by the odd shell at \(n/q^2\), giving
\(wS/q^2+(1-w)\Delta/q\) after division by four.

At \(p\), the bad block norms lie on the axes in the radical quotient.
Their indicator is the principal base point and has mean
\(\nu=2w/(p+1)\). Its projections onto the constant and quadratic
norm characters each have coefficient \(\nu\).
The bad count divided by four is at most
\((1+(n/p))\nu S+(1-2\nu)\Delta\), by the same self-correlation.
[Local normalizer and radical correlations](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md#small-bad-prime-counts-and-the-ramified-filter)

Let \(M(X)\) be the maximum over \(1\le v\le X\) of the divisor product
\(\prod_{q\mid v,\ (q/p)=1}(v_q(v)+1)\). By class number one,
\(r_{Q_p}(v)\le2M(X)\) in that range. For an odd inert prime \(q>P\),
a bad block is \(q\) times a nonzero block of norm at most \(n/q^2\).
Since \(Q_p\ge Q_{11}\), the lattice bound (M11.17) gives at most
\(12n/(5q^2)\) choices. Counting the complementary block and the two
orientations bounds the bad count divided by four by
\(12M(n)n/(5q^2)\).
[Binary lattice bound (M11.17)](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md#closing-the-finite-core-sharp-bounds-and-a-complete-certificate)

Write \(B=\sum_q\mu_q\), \(\Pi(P)=\#\{q:q\le P\}\) and
\(T(P)=\sum_{q>P}q^{-2}\), with all sums over odd inert primes.
Using \(S\ge3n/2\), subtraction of the bad counts proves positivity once
\[
 B+\frac{8}{5w}M(n)T(P)+
   \frac{2-w-2\nu+\Pi(P)}{w}E < 1-\frac4{p+1}.                  \tag{1}
\]
For selected odd inert \(q\le P\) dividing \(n\), the left side can
replace the contribution \((\mu_q,1)\) to the mass and error
coefficient by \((0,0)\) at valuation one and \((q^{-2},q^{-1})\)
at larger valuations. Every surviving point has two positive core
norms. Split factors in its coordinates are permitted because it is
the numerical values that must belong to \(C_p\).

## Infinite estimates and completion at −19

For the four primes, the rational checks give \(M(n)\le7n^{1/6}\),
\(E\le16n^{-1/3}\), and \(B<27/50\).
For the first two bounds, the sixth-power Euler factors to maximize at
prime \(q\) and exponent \(a\) are
\[
 (a+1)^6/q^a,\qquad
 (a+1)^6 q^{5a}/\sigma(q^a)^6,
\]
respectively. The second product also has the dyadic factor bounded
by \((3/2)^6\), since \((t+3)2^{-t/6}/3\le3/2\).
The primes \(q\ge64\) contribute at most one by \(a+1\le2^a\).
For smaller primes, \(0\le a\le10\) suffices since the ratio of
successive factors for \(a\ge10\) is bounded by
\((12/11)^6/q<1\) for odd \(q\ge3\).
For \(B\), the certificate sums the odd inert primes through
\(R=65536\) with upward rounding, then bounds the tail by summing over
the \((p-1)/2\) odd inert progressions modulo \(2p\).

For general \(P>0\) those progressions give
\(\Pi(P)\le P/4+81\) and \(T(P)\le1/(4P)+81/P^2\).
Take \(P=n^{1/4}\) and use \(w\ge1/27\). At \(n\ge2^{120}\) the
left side of (1) is bounded by
\[
 \frac{27}{50}+\frac{216}{5}\,7
       \left(\frac1{4\cdot2^{10}}+\frac{81}{2^{40}}\right)
   +27\cdot16\left(\frac1{4\cdot2^{10}}+\frac{83}{2^{40}}\right)
   <\frac45\le 1-\frac4{p+1}.
\]
This tail applies to every allowed target, including ones with split
prime factors.

At \(p=19\), it suffices to handle two restricted families above the
scan bound. One family consists of even targets supported on inert
primes. The other is \(n=4q\) or \(20q\) for a single split prime \(q\).
In the first family the finite Euler products bound \(En^{1/3}\)
by constants indexed by the valuations of 3 and 13, grouped as
zero, one or at least two. Every such constant is less than 5.
For the second family \(q>2^{12}\), and the factor at \(q\) is at
most \(2q^{-1/6}\le1/2\) before taking the sixth power.
The factor at 5, when present, is at most 2. Thus
\(En^{1/3}\le3/2\) on the second family.

For \(2^k\le n<2^{k+1}\), \(30\le k\le71\), bound \(E\) using the
corresponding constant divided by \(\lfloor\sqrt[3]{2^k}\rfloor\).
The divisor envelope at \(2^{k+1}-1\) comes from the initial split
prime prefixes with nonincreasing exponents. Reassigning decreasing
exponents to increasing primes weakly decreases the integer and
preserves its divisor product, proving completeness of the envelope.
The checker certifies (1) by choosing a prime cutoff in every one of
the 420 bin and profile pairs. It includes the valuation corrections
for 3 and 13 and adds the residue-progression tail to each rounded
finite tail sum. For \(n\ge 2^{72}\) in these two families use
\(E\le5n^{-1/3}\), \(M(n)\le7n^{1/6}\), \(P=n^{1/4}\) and the nine
odd inert progressions modulo 38. Their sharper bounds
\(\Pi(P)\le9P/38+9\) and \(T(P)\le9/(38P)+9/P^2\) verify (1) on the
whole remaining range by monotonicity of the negative powers of \(n\).
[Rational checker](verify_heegner_core_bounds.py)

If a core exception outside \(\{4,20\}\) exists, take the least one.
It exceeds \(2^{30}\) by the scan. An inert-supported target is covered
above. Otherwise divide by a split prime \(q\). The quotient remains
even and a unit at 19. If represented, it lifts to the target.
If unrepresented, minimality forces the quotient to be 4 or 20, and
those prime multiples were covered by the second family. This
exhausts the −19 classification.


## Sharper tails at −43, −67 and −163

For these discriminants count only odd blocks in the tail. Let \(H(X)\)
count all nonzero \(Q_p\)-vectors of norm at most \(X\).
For integral \(y\) with \(|y|\le2\sqrt{X/p}\), the row length is
\(2\sqrt{X-py^2/4}\), and the integral point count differs from that
length by at most one. Comparing the decreasing half-row sum with
its integral shows
\[
 A X-b_0\sqrt X-2\le H(X)\le A X+b_0\sqrt X,\qquad
 A=2\pi/\sqrt p,\quad b_0=2+4/\sqrt p .
\]
The norm is even exactly on the even coordinate vectors, so the odd
count is \(H(X)-H(X/4)\). It is at most \(aX+b\sqrt X+2\), where
\(a\ge3\pi/(2\sqrt p)\) and \(b\ge3+6/\sqrt p\).
The following rational choices follow from \(\pi<22/7\).

| \(p\) | \(a\) | \(b\) | First inert bin \(k\) | First \(k\) of the infinite bound |
| --- | ---: | ---: | ---: | ---: |
| 43 | \(18/25\) | 4 | 32 | 72 |
| 67 | \(29/50\) | 4 | 35 | 84 |
| 163 | \(37/100\) | \(7/2\) | 43 | 84 |

Put \(F=\prod_{q} \sum_{j=0}^{v_q}q^{-j}\), over the tracked odd inert
primes. Use the four smallest at \(p=43,67\) and the six smallest at \(p=163\).
Here \(v_q\) is the truncated valuation, at most two.
Then \(S\ge(3/2)nF\). For the other profile take \(F=1\).
In (1) replace the large-prime term by the sharper bound
\[
 \frac{2M(n)}{3wF}
 \left(aT(P)+\frac{b}{\sqrt n}U(P,\sqrt n)
                    +\frac{2\Pi(\sqrt n)}n\right),
 \qquad U(P,Y)=\sum_{P<q\le Y}q^{-1}
\]
over odd inert primes. This follows by replacing the count of the
bad block divided by its inert prime with the odd lattice count.

For the inert-supported targets the certificate uses \(3^4\) valuation
profiles at \(p=43,67\) and \(3^6\) at \(p=163\). It adjusts the main
mass and error terms at target-dividing primes as in (1).
For each profile it bounds \(En^{1/3}\) by the corresponding finite Euler
maximum. The other profile treats \(n=eq\) with \(e\) one of the
observed gaps and \(q\) a split prime. It starts at the scan bound
in each case, including \(2^{34}\) at \(p=163\).
Here \(q>e\), and if \(m_0\) is the odd part of \(e\), then
\[
 En^{1/3}\le\frac{3}{2}\frac{\tau(m_0)m_0^{5/6}}{\sigma(m_0)}\,2q^{-1/6}.
\]
The checker maximizes the sixth power over odd parts up to the
largest candidate divided by two, using the scan bound divided by
the largest candidate to bound \(q\) from below.

The checker uses the same exact divisor frontier and square- or
cube-root floors to bound the negative powers of \(n\) in each bin.
For \(U(P,\sqrt n)\) it sums through \(R\) with upward rounding.
At larger primes, summing \(1/q\) over each odd inert progression
from \(R\) to \(2^{\lceil(k+1)/2\rceil}\) and using \(\log 2<7/10\)
bounds the remaining term. It checks 3280, 4018 and 29939 bin and
profile inequalities for \(p=43,67,163\), respectively.

For the infinite estimates, let \(c=(p-1)/2\), \(r=c/(2p)\), and
\(x=n^{1/12}\). Take \(P=x^3/2\). Progression bounds give
\[
 T(P)\le 2r/x^3+4c/x^6,\qquad
 U(P,\sqrt n)\le 2c/x^3+r(1+3x),\qquad \Pi(Y)\le rY+c .
\]
The middle bound uses \(\log(2x^3)\le1+3x\).
The Euler certificates give \(En^{1/3}\le7,9,14\) and
\(M(n)n^{-1/6}\le4,3,2\) on the respective reduced families.
After substitution every \(x\)-dependent term decreases.
Checking at the first \(k\) of the infinite bound covers all
larger reduced targets.

The scans for −43 and −67 meet the first bins in the table.
The same least-exception argument as for −19 completes both
classifications. At −163, a least missing target outside the
retained exception list could only be inert-supported and lie
strictly between \(2^{34}\) and \(2^{43}\). Split-prime multiples
of the listed gaps are covered from the scan bound onward.
The next section certifies the remaining inert-supported range.


## Completing −163 without scanning the whole interval

The sharper tail and the old finite scan leave even inert-supported
targets \(2^{34}\le n<2^{43}\). The left endpoint can be included
again. Write \(n=2^t m\) with \(m\) odd. On this domain the normalized
error is exactly
\[
 E=\frac{(t+3)D(m)}{3\sqrt n},
 \qquad D(m)=\frac{\tau(m)m}{\sigma(m)} .
\]
In addition to the six valuations already used for the tail, distinguish
the sign \(\chi=(n/163)=(-1)^{t+\Omega(m)}\).
When \(\chi=-1\) the bad ramified axes are empty, since one nonzero
square on an axis would make the target a square. Both their main
term and their spectral error can therefore be omitted.

Take a dyadic bin \([L,U]=[2^k,2^{k+1}-1]\).
For each of the \(3^6\) profiles, write \(b_v\) for the adjusted upper
bound for \(B\), \(F_v\) for the product used to bound \(S/n\), and
\(A_{v,\chi}\) for the adjusted constant in the spectral coefficient.
Without any corrections \(A_{v,\chi}=1-w\) plus \(1-2\nu\) for
\(\chi=1\). Each target-dividing tracked prime replaces its term
\(1\) in the coefficient by \(0\) or \(1/q\), as in (1).
If \(P\) is the \(i\)-th odd inert prime and includes the six tracked
primes, the sharpened sieve therefore certifies representation once
\[
 b_v+\frac{V_{k,i}}{F_v}+(A_{v,\chi}+i)\frac{E}{w}
   < R_\chi,\qquad R_1=1-\frac4{164},\quad R_{-1}=1.
\]
Here \(V_{k,i}\) is the upper bound for the odd lattice tail from the
preceding section with \(F=1\), using the exact divisor envelope at \(U\)
and downward square-root rounding at \(L\). In its count term use
\(2\Pi(\sqrt n)/n\le 1/(2\sqrt n)+162/n\).
The same finite prime sum and progression bounds give the two other
tail terms.

The checker chooses an \(i\) by exact rational comparisons and retains
the positive integer threshold
\[
 K_{k,v,\chi}
  =\left\lfloor 2^{20}\,3\lfloor\sqrt L\rfloor\, w
       \frac{R_\chi-b_v-V_{k,i}/F_v}{A_{v,\chi}+i}\right\rfloor .
\]
Every possible failure satisfies
\(2^{20}(t+3)D(m)\ge K_{k,v,\chi}\).
Only such targets need an individual witness. Any selected \(i\)
giving a positive threshold is sufficient for this implication.

To see that the recursion needs primes only through \(2^{25}\),
suppose a remaining target has an odd prime \(q>2^{25}\).
Then \(q\) occurs once, since \(q^2>U\). The odd integer \(u=m/q\)
is inert-supported and at most \(\lfloor U/(2^t(2^{25}+1))\rfloor\).
Since \(D(m)\le\tau(m)=2\tau(u)\), the exact inert divisor envelope
on that small interval rules out such \(q\) whenever the threshold
inequality is required. The checker verifies this against the minimum
of the profile thresholds for every \(k,t\).

Enumerate the six small-prime exponents first, preserving the base-three
code of the truncated valuations, then enumerate larger inert primes
in increasing order. For a partial odd product \(u\), the needed
multiplier of \(D(u)\) is bounded below by using the smaller threshold
over the two possible signs. If its ceiling is \(h\), any remaining
factor \(v\) must have \(\tau(v)\ge h\). Hence it has at least
\(\lceil\log_2h\rceil\) prime factors with multiplicity, by
\(D(v)\le\tau(v)\le2^{\Omega(v)}\). If the next prime is \(q\), their
product is at least \(q^{\lceil\log_2h\rceil}\). The recursion stops
the branch if this already exceeds \(U/(2^t u)\).
It uses at least one prime factor whenever it considers a nonempty
continuation. Every unpruned terminal product is checked against
the threshold for its actual sign. Unique factorization makes these
branches disjoint and exhaustive.

The resulting shortlist has \(15\,406\,572\) targets.
For each, the native checker finds a positive core summand \(a\) from
an exact small sieve. It strips recognized split prime factors from
\(n-a\). The residue character remains one, and the remaining factor is
either one or a certified prime, so all its prime factors split as well.
Primality is obtained by Eratosthenes in the retained range and seven
strong prime-base tests above it. Jaeschke's exact bound
\(341550071728321>2^{43}\) for these bases makes the latter test
deterministic on this domain.
[Jaeschke, Section 3, case SC2, printed pages 924–925](sources/jaeschke-strong-pseudoprimes-1993.pdf)

The largest shorter summand used was \(5351\).
The top bin \(k=42\) had just one shortlisted target,
\(8236350794672=499+8236350794173\), with the second summand
\(47\cdot71\cdot2468190229\), a product of split primes.
An independent segmented factorization of two 100000-integer sample
windows recovered the same 52 and 9 shortlisted targets.
Trial division independently verified their witnesses and the top-bin
witness. These checks are additional internal controls. The exact
completion follows from the rational pruning and complete traversal.
[Completion receipt](heegner_core_completion_checks.json)


## A common effective bound, including −7

For \(p=7\) the same quaternionic construction, unit group,
mass and radical quotient apply. The only change in the initial
shell is at 2. Since 2 splits in the quadratic field, it is allowed in
both summands, and the operator is \(T_n\) itself for every
\(n\ge2\) prime to 7. Thus \(w=1\), \(S=\sigma(n)\),
\(\Delta=\tau(n)\sqrt n\), and \(\nu=1/4\).
The same odd inert normalizer and ramified-axis correlations now
give (1) with \(8/(5w)\) replaced by \(6/w=6\).
For that large-prime term, count only points surviving the
ramified filter, so neither block is zero. The bound of \(6X\)
nonzero binary lattice vectors of norm at most \(X\) follows by
comparison with the Eisenstein form, as in (M11.12).
Use \(S\ge n\) in this branch. No restriction on the parity of \(n\)
remains.

At \(p=7\) the certificate gives \(B<47/100\),
\(M(n)\le8n^{1/6}\), and \(E\le20n^{-1/3}\).
The last bound uses the exact \(\sigma(n)\) denominator and the same
sixth-power Euler factors as at the other primes, now over all
primes including 2. Exponents through 10 suffice also at 2 since
\((12/11)^6<2\). The three odd inert progressions modulo 14 give
\(\Pi(P)\le 3P/14+3\) and \(T(P)\le3/(14P)+3/P^2\).
With \(P=n^{1/4}\), for \(n\ge2^{120}\) the left side of (1) in
this branch is at most
\[
 \frac{47}{100}
 +48\left(\frac3{14\cdot 2^{10}}+\frac3{2^{40}}\right)
 +20\left(\frac3{14\cdot 2^{10}}+\frac7{2\cdot 2^{40}}\right)
 <\frac12.
\]
This is the required strict inequality at \(p=7\).
For the other four primes the previously derived bound gives the
same cutoff. Thus the five-case effective theorem uses one shell
and sieve with the local choice of whether to filter at 2.
The exact classification for \(p=7\) uses Theorem M7.1 of the
published argument. The other four exact classifications follow
from the certified finite ranges and split-prime descent here.

## Reproduction

From the repository root, use a C++20 compiler and Python 3 with
assertions active.

```sh
g++ -std=c++20 -O3 -march=native -Wall -Wextra -pedantic \
  research/arithmetic/verify_heegner_core_sumsets.cpp -o /tmp/verify_heegner_core_sumsets
/tmp/verify_heegner_core_sumsets 1073741824 19
/tmp/verify_heegner_core_sumsets 4294967296 43
/tmp/verify_heegner_core_sumsets 34359738368 67
/tmp/verify_heegner_core_sumsets 17179869184 163
python3 research/arithmetic/verify_heegner_core_bounds.py
g++ -std=c++20 -O3 -march=native -fopenmp -Wall -Wextra \
  research/arithmetic/verify_minus163_core_bridge.cpp -o /tmp/verify_minus163_core_bridge
OMP_NUM_THREADS=8 python3 research/arithmetic/verify_heegner_core_completion.py \
  --native /tmp/verify_minus163_core_bridge
```

The native scan represents odd integers by packed bits. It deletes
multiples of every non-split prime, starting from an exact prime sieve.
Every shift in its union is by an actual core summand at most 5000.
For each admissible target still uncovered, it tests every possible
smaller odd summand up to half the target. For each resulting gap it
tests every possible mixed pair. The checker cross-checks the sieve
against primitive-coordinate enumeration through 10000.

The retained execution used GCC 13.3.0 and Python 3.12.3.
A separately implemented internal coordinate-and-sieve scan agreed
on all exception lists through \(2^{25}\). The split-core dilation
and mixed tests also reproduce the previously retained primitive
lists through \(10^6\). The finite scan and rational checker have
different scopes from the written spectral argument.
[Exact receipt](heegner_core_sumset_checks.json) ·
[Previous primitive scan](primitive_norm_exact_sumset_checks.json)
