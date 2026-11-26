# Short projective intervals with small radical

The short-interval radical bound extends to reduced positive
rationals. The analogue of a large full divisor is now a divisor
of a cross determinant. For primitive abc triples this controls
entropy inside very short intervals of a pair ratio and gives
a real localization alternative to missing valuation rank.

## A projective short-interval bound

Fix \(A,H,h_0,\eta>0\). There are fixed \(\theta,L_0>0\)
with the following property. Let \(0<L<L_0\) be fixed
before \(X\to\infty\), and let \(h_1,h_2\in[0,H]\) be
height caps with \(h_1+h_2\ge h_0\). For every
\(\varepsilon>0\), uniformly in \(h_1,h_2\) and \(K>0\),
\[
 \#\left\{u/v:
 \begin{array}{l}
 u,v\in\mathbb Z_{>0},\ (u,v)=1,\quad
 u\le X^{h_1},\quad v\le X^{h_2},\\
 \operatorname{rad}(uv)\le X^{AL},\quad
 \lvert (u/v)/K-1\rvert\le X^{-(h_1+h_2)+\theta}
 \end{array}\right\}
                   \ll_{A,H,h_0,\eta,L,\varepsilon}
                           X^{\eta L+\varepsilon}.              \tag{PI1}
\]
The same statement holds with a fixed constant multiplying the
strip width, by reducing \(\theta\). The height caps may vary
with the shell. In particular two height caps of one give
relative width \(X^{-2+\theta}\).

Put \(b=\min(A/4,\eta/16)>0\). In a candidate \(u/v\), there
are at most \(M=\lfloor A/b\rfloor\) distinct primes
\(p_i\ge X^{bL}\). Write
\[
                          u/v=V\prod_{i=1}^m p_i^{a_i}
\]
in increasing order of these primes, using signed exact powers
and a reduced positive cofactor \(V\) supported on the smaller
primes. All \(p_i\le X^{AL}\) and the nonzero exponents
satisfy \(|a_i|\le H/(bL)\). At every fixed budget there
are finitely many possibilities for these signed order lists,
uniformly in the two height caps. The displayed presentation
is unique for each reduced fraction.

Choose \(s>0\) so that \((M^2+M)s<\eta/4\). Fix \(E>0\) that
satisfies the strictly wider-height margin in (RUW1) for
every \(1\le m\le M\), with narrow cap parameters \(A+1\),
strip precision \(\tau=h_0/2\), and the indicated \(s\).
Also take \(E<h_0/8\). Take \(L_0\) below the budget
thresholds in these finitely many applications. Choose
\[
                    0<\theta<\min\{h_0/4,E/(2(M+1))\}.
\]
For \(H(V)\le X^E\) and \(m\ge1\), (RUW1) applies in
dimension \(m\) with the fixed order list and with exact
monomial fibre size one. Its cost fits \(X^{\eta L/4+\varepsilon}\).
If \(m=0\) and \(H(V)\le X^E\), two distinct such rationals
would have a nonzero integer cross determinant, yet the
relative strip bounds its magnitude by
\(O(X^{-(h_1+h_2)+\theta+2E})<1\).
[Real strip estimate (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

If \(H(V)>X^E\), choose the numerator or denominator of \(V\)
whose height exceeds \(X^E\). On that side group complete
prime powers greedily until the next prime would push the
radical of a group above \(X^{2bL}\). Each completed group
other than the last has radical larger than \(X^{bL}\),
since every prime on the chosen side is smaller than
\(X^{bL}\). There are at most \(M+1\) groups. One yields a
full divisor \(F\ge X^{E/(M+1)}\) whose radical is at most
\(X^{2bL}\). Record its value and its side. There are
\(O(X^{2bL+\varepsilon})\) such records by enumerating the
radical and then the possible powers at height at most
\(X^H\).

For two candidates with the same record,
\(F\mid u_1v_2-u_2v_1\). The relative strip gives
\[
        |u_1v_2-u_2v_1|
            \ll X^{-(h_1+h_2)+\theta} u_2v_1
            \le O(X^\theta).
\]
For sufficiently large \(X\), this is smaller than \(F\),
so the two reduced ratios coincide. Thus each record costs
at most one candidate. Since \(2b\le\eta/8\), this proves
(PI1) with arbitrary inner exponent allowances.
[Integer short-interval analogue (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)

## Real localization in an approximate valuation packing

Work on a dyadic shell of positive primitive abc triples with
total radical cap \(X^L\). Write \(Z=U/V\) for a prescribed
pair ratio. First use dyadic complete-height bins with
caps \(X^{h_1},X^{h_2}\) on \(U,V\). At least one of
\(U,V\) is at least \(X/4\), since two of the three letters
are at that height. One can use \(H=1\), \(A=1\), and
\(h_0=1/2\) in (PI1) for large \(X\).

Within each fixed height bin, divide \(\log Z\) into equal
intervals of length \(\tfrac14 X^{-(h_1+h_2)+\theta}\).
Let \(J\) record the interval and the complete-height bin.
By (PI1), each value of \(J\) admits at most
\(O(X^{\eta L+\varepsilon})\) possible exact ratios. Each
distinct ratio recovers a primitive row on the prescribed
roles. For every prime projection \(Y_S=(v_p(Z))_{p\in S}\),
\[
                  H(Z\mid Y_S)
                  \le H(J\mid Y_S)
                       +(\eta L+\varepsilon)\log X+O(1).
                                                               \tag{PI2}
\]
In particular a selector can localize the ratio in a small
dictionary of numerical bins without determining its
unobserved valuation coordinates. If a separate valuation
kernel estimate is available, one may use the smaller of
(PI2) and (AV3) for that selector.
[Entropy and kernel-rank bounds (AV1–AV3)](ABC_APPROXIMATE_VALUATION_SEPARATORS.md)

For a fractional packing of weight \(P\), suppose its
weighted mean \(H(J\mid Y_S)\) is at most
\(\omega L\log X\). Substitution in (AV1) and the ratio
case of (AV2) yields, in the balanced limit with individual
radical caps \(1/3+o(1)\) in base \(X^L\),
\[
                     d\le \frac{2}{3P}+\omega+\eta+o(1),
\]
with inner errors taken before the budget limit. Any fixed
\(P>2\) and \(\omega+\eta<1/3-2/(3P)\) therefore give a
strict saving. One sufficient condition is a uniform bound
of \(X^{\omega L}\) on the number of projective bins allowed
after each of the projected prime records. The same
numerical argument bounds the conditional row entropy when
the packings are available in separate recorded fibres. The
entropy of the conditioning record is then added for the
unconditional count. The extraction of such localizing projections from a
general balanced family is a separate step.
