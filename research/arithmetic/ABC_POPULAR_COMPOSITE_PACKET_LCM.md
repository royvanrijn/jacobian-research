# Exponent-popular composite packet LCMs

In a fixed original fibre, a bounded list of popular congruence
classes cannot resolve almost the full numerical height of a moving
complete letter on many rows. The divisors defining the classes can
vary with the row, overlap, and contain arbitrarily small primes.
The same restriction holds in an actual common-target completion
star of several fixed complete summands.

## A bound for row-dependent popular divisors

Fix \(C_0,B,r>0\) and a positive integer \(K\). Let \(n=X^L\)
and work with positive primitive rows \(A+B'=C\) in
\(X/2<C\le X\), whose three individual radicals are at most
\(n^{C_0}\). Fix one complete letter and a different moving
letter \(W\). Write \(\mathcal F\) for any set of \(m\) rows
of this fibre. The complete value \(W\) identifies its row.
Assign each such value a nonzero integer \(R(W)\) with
\(|R(W)|\le X^B\). Require the transfer property
\[
 D\mid R(W),\quad W'\equiv W\pmod D,\quad W,W'\in\mathcal F
       \ \Longrightarrow\ D\mid R(W').
                                                               \tag{PL1}
\]
In (PL1), membership means that the corresponding rows belong to
\(\mathcal F\). An integer polynomial \(R(W)\) satisfies (PL1)
regardless of how its coefficients depend on the prescribed fibre.

For a tested value \(W\), allow any list of \(k\le K\) positive
divisors \(Q_i\mid R(W)\). The list is eligible if
\[
 \begin{split}
  \operatorname{lcm}_i Q_i&\ge X^{1-\theta},\\
  \prod_{i=1}^k
   \frac{\left|\{W'\in\mathcal F:W'\equiv W\pmod {Q_i}\}\right|}{m}
                    &\ge n^{-p}.
 \end{split}                                                    \tag{PL2}
\]
In particular the second condition holds when the \(i\)-th
class has relative frequency at least \(n^{-\tau_i}\)
and \(\sum_i\tau_i\le p\). Relative frequencies always refer
to the prescribed \(\mathcal F\).

**PL3.** There are \(\theta,\sigma,L_0>0\) depending
only on \(C_0,r\) such that for every fixed \(0<L<L_0\)
and every \(p\ge0\), the number of values \(W\) admitting
(PL2) is at most
\[
             O_{C_0,B,K,r,L,\varepsilon}
                      \bigl(n^{r-\sigma+p}X^\varepsilon\bigr).
                                                               \tag{PL3}
\]
In particular, if \(m\ge n^r\), these values form a
power-small fraction when \(p<\sigma\) with a fixed margin.
The exponent threshold is independent of the size of the
primes dividing the packets. All heights here tend to
infinity with \(L\) fixed.

To prove this, use the rank-one small-feature estimate
at exponent \(r\), and denote its feature height exponent and
count gap by \(\theta_0,\sigma\). Choose
\(0<\theta<\min(\theta_0,1)/2\). For an eligible \(W\) choose
one list independently of any witness sampling. This can be
done deterministically, by taking the first eligible list in
a fixed order. Pad it to length \(K\) with \(Q_i=1\).
Sample \(K\) independent uniform anchor rows
\(W_1,\ldots,W_K\) from \(\mathcal F\), with replacement.
The probability that every anchor lies in the corresponding
class in (PL2) is at least \(n^{-p}\).

Now fix the sampled anchors. On that event every \(Q_i\)
divides \(R(W_i)\) by (PL1). Enumerate these \(K\)
divisors at \(O_{B,K,\varepsilon}(X^\varepsilon)\) cost
uniformly. The \(K\) congruences \(W\equiv W_i\bmod Q_i\)
are either inconsistent or specify a single class modulo
\(D=\operatorname{lcm}_iQ_i\). For \(D\ge X^{1-\theta}\),
the values in that class between \(1\) and \(X\) are
\(W=W_0+D z\) for a fixed integer \(W_0\) and integer
\(z\) of absolute value at most \(1+X^\theta\), by taking
\(W_0\) to be a least nonnegative representative. The rank-one
estimate counts \(O_{C_0,r,L}(n^{r-\sigma})\) such rows.
This uniform bound is valid after conditioning on all the
anchors and divisors. Averaging over the anchors and using the
lower success probability for each tested value proves (PL3).

The argument allows repeated anchors and even an anchor equal to
the tested row. Such coincidences do not affect divisibility or
the rank-one bound. If witnesses must instead be distinct,
different from the tested row, and avoid any previously fixed
list of bounded size, the same count follows under \(m\ge n^r\) and
\(p<r\) with a fixed margin. Every factor in the product in
(PL2) is then at least \(n^{-p}\). Each class has
\(n^{r-p}\) or more members, so deleting the forbidden points
and the previously chosen witnesses lowers its share by at most
a factor two for large \(X\). The success probability remains
at least \(2^{-K}n^{-p}\). The packet list, and its required
frequencies or modes, may depend on the fixed fibre and tested
row, independently of the sampled witnesses.
[Rank-one fixed-fibre estimate](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

## The three original letters and completion stars

For a fixed complete target \(C\), take \(W=A\) and
\(R(W)=W(C-W)C\). An eligible packet may combine prime
powers from the two actual summands and the target.
For a completion star of fixed positive summands
\(U_1,\ldots,U_s\), use an arbitrary set of complete targets
\(W\) for which every \((U_i,W-U_i,W)\) is an original
primitive row with all three caps. Set
\[
                       R(W)=W\prod_{i=1}^s(W-U_i).
\]
Its height is at most \(X^{s+1}\). The rows attached to
\(U_1\) supply the fixed fibre for (PL3). Thus the theorem
restricts how much depth can be recovered from popular
divisor classes of all the actual moving letters in a star.
It applies with \(s\) fixed and does not require the
factors of the displayed product to be coprime.

There is also a depth statement with constant popularity.
Fix \(0<\alpha\le1\). For a row \(W\) and a prime power
\(q\mid R(W)\), call the depth \(q\) popular if the
relative frequency of \(W'\equiv W\bmod q\) in \(\mathcal F\)
is at least \(\alpha\). Form the product \(D_\alpha(W)\)
of the largest popular powers of each prime. The set of
values with \(D_\alpha(W)\ge X^{1-\theta/2}\) satisfies
\[
           O_{C_0,B,r,\alpha,L,\varepsilon}
                         \bigl(n^{r-\sigma}X^\varepsilon\bigr).
                                                               \tag{PL4}
\]
For the proof sample \(J\) anchors with replacement, where
\(J\) depends on \(B,\theta,\alpha\), and choose for each anchor
the divisor \(\gcd(|W-W_i|,|R(W_i)|)\). This includes the
usual convention \(\gcd(0,a)=a\). A popular depth fails to
divide the least common multiple of the chosen divisors only
if all the anchors miss its residue class, by (PL1). The
expected missing logarithmic mass of \(D_\alpha(W)\),
in base \(X\), is at most \(B(1-\alpha)^J\).
Choose \(J\) so this is at most \(\theta/8\). By Markov's
inequality, with probability at least \(3/4\) the missing mass
is at most \(\theta/2\). The resulting least common multiple
is then at least \(X^{1-\theta}\). For each fixed anchor
tuple, enumerate the \(J\) divisors and apply the same
rank-one bound as above. Averaging and dividing by the
constant success probability gives (PL4).

For example, take \(\alpha=n^{-\tau}\) but restrict
to prime powers \(q\) whose prime is at least \(n^\beta\)
and belongs to one of the original moving letters in a
fixed completion star. There are at most
\((s+1)C_0/\beta\) such primes on each tested row.
Listing their maximal popular powers as separate packets
therefore puts this exponent-scale version directly under
(PL3) for fixed \(\beta>0\). The bounded composite
packet form (PL3) also allows many small primes per packet.

These estimates exclude near-total numerical LCMs of popular
absolute congruence classes. A target-fibre tail argument
still needs to force both the popularity and the numerical
LCM at the requisite scales, or to bound the complementary
families. In particular, (PL3) does not imply that the physical
shifts have low signed radical cost as in (SCM8).
[Signed costs and actual target codegrees](ABC_COMPRESSED_SHIFT_TARGET_TAILS.md#endpoint-support-and-original-target-discounts)


## Forced high-valuation tails in common completion stars

Fix \(C_0,r>0\) and use \(\theta,\sigma,L_0\) from (PL3).
Let \(J\ge1\) be fixed, and suppose
\(F_1(W),\ldots,F_J(W)\) are prescribed integral affine
functions of slope \(1\) or \(-1\) with positive values
at most \(X\) on the fixed original fibre. For \(R>0\)
let \(Q_j(W)\) be the product of the full prime powers
\(p^e\parallel F_j(W)\) with \(Le\ge R\).

**PL5.** With \(L\) sufficiently small and fixed as in
(PL3), the number of \(W\) with
\(\operatorname{lcm}_j Q_j(W)\ge X^{1-\theta}\) is
\[
     O_{C_0,r,J,R,L,\varepsilon}
                    \bigl(n^{r-\sigma+2/R}X^\varepsilon\bigr).
                                                               \tag{PL5}
\]
The exponents \(\theta,\sigma\) and the required choice of
cutoff \(R\) for a prescribed saving are independent of the
fixed number \(J\) of affine forms. The assertion allows the
tail factors of the forms to overlap.

Indeed, take full prime powers of the least common multiple
until the product \(D\) first reaches \(X^{1-\theta}\). A
single chosen power is at most \(X\), so \(D\le X^2\).
Its radical is at most \(n^{2/R}\) by the depth condition.
Record the radical, its powers, and for each selected prime
one form attaining the selected depth. There are at most
\(n^{2/R}X^\varepsilon\) records, since the powers and the
allocations into \(J\) forms cost a subpower factor at fixed
\(J\). Each affine equation \(F_j(W)\equiv0\bmod p^e\)
fixes \(W\bmod p^e\). The recorded prime powers have
distinct primes, so they specify a single progression
modulo \(D\) when consistent. Apply the rank-one bound in
that progression. It is not necessary to enumerate all
the high tails of all the forms: stopping at one height's
worth of numerical modulus keeps the charge \(2/R\)
independent of \(J\).

This gives a dichotomy in an actual summand–target
projection. Let \(\mathcal E\) be a set of original
primitive rows with the same height shell and individual
caps. Let its source dictionary on \(A\) have size at most
\(n^a\), and fix \(u,\omega>0\) with \(\omega\le1\).
For each complete target \(C\), write \(\mathcal A_C\) for
its \(A\)-neighbours and \(m_C=|\mathcal A_C|\). For a
neighbour \(x\) let \(Q_C(x)\) be the full high tail of the
actual opposite letter \(C-x\), at cutoff \(Le\ge R\) as
above. Define the mean numerical tail weight and the mean
pairwise common weight by
\[
       \omega_C=\frac1{m_C}\sum_{x\in\mathcal A_C}
                                     \log_X Q_C(x),\qquad
       \xi_C=\frac1{m_C^2}
                 \sum_{x,y\in\mathcal A_C}
                               \log_X\gcd(Q_C(x),Q_C(y)).
\]
The second sum includes the diagonal.

**PL6.** Let \(s=\lceil 3/\omega\rceil\). The number of
complete targets in the cohort
\[
       m_C\ge n^u,\qquad \omega_C\ge\omega,\qquad
       \xi_C\le\frac{\omega}{s-1}
\]
is at most
\[
       O_{C_0,r,s,R,u,L,\varepsilon}
            \bigl(n^{s(a-u)+r-\sigma+2/R}X^\varepsilon\bigr).
                                                               \tag{PL6}
\]
In particular, for any fixed \(\eta>0\), one can choose
\(r>0\), a sufficiently large fixed \(R\) and then a
small-budget threshold so that the exponent is at most
\(s(a-u)+\eta\), by choosing the inner allowance at
each fixed budget. The choices of \(r,R\) and the budget threshold
are independent of the fixed star size \(s\), apart from
the inner counting constants.

For each target in this cohort sample \(s\) independent
uniform neighbours. The numerical weight of the least
common multiple of their tails is at most \(s\), and its
expectation is at least
\[
               s\omega_C-\binom{s}{2}\xi_C
                         \ \ge\ \frac{s\omega}{2}\ \ge\ \frac32 .
\]
The first inequality holds separately at every prime
power depth by the first two terms of inclusion–exclusion.
A constant fraction, depending only on \(s\), of the
sampled stars therefore have numerical LCM weight at least
one. Removing repeated neighbours loses at most
\(\binom{s}{2}/m_C\) of the sample. Thus each such target
contributes \(\gg_s m_C^s\) ordered stars with distinct
sources and large LCM.

Now fix the ordered source tuple \((x_1,\ldots,x_s)\).
Any actual common target belongs to an original
fixed-letter fibre at \(x_1\) and gives the affine values
\(W-x_1,\ldots,W-x_s\) with \(W=C\). Apply (PL5) to count
the stars of large LCM for that tuple. There are at most
\(n^{sa}\) tuples. Summing and dividing by \(n^{su}\)
proves (PL6).

For a target cohort of exponent \(c\) and filtered edge
exponent \(d=c+u\), the comparison is with
\(s(a+c-d)+\eta\). A fixed margin \(c>s(a+c-d)+\eta\)
gives a power saving for the low-common-weight branch.
The other branch, with \(\xi_C>\omega/(s-1)\), has a
constant fraction of distinct pairs whose common high tail
has numerical weight at least \(\omega/(2(s-1))\).
For each such pair let \(h=|x-y|\). This common high tail
divides \(h\), has radical at most \(n^{1/R}\), and is
coprime to the normalization factor \(\gcd(x,h)\) of the
\(A\)-side difference by original primitivity. It therefore
survives on that normalized shift. Taking the complete powers
of the same primes in the normalized shift preserves the radical
and only increases their depths and numerical weight. This gives
shifts with a prescribed minimum depth and numerical weight
of a full divisor, including shifts of large raw radical.
The conclusion requires mean numerical mass in the high
tails. It gives no upper bound on the signed shift cost of
the remaining primes and no saving for general balanced
profiles.


## A sharp occupancy criterion retaining the target tail

Keep the notation of (PL6), let \(s\ge1\) be any fixed
integer, and use \(\theta,\sigma\) from (PL5).
Let \(P_C\) be the full high tail of the actual target \(C\)
at the same cutoff \(Le\ge R\), and write
\(\gamma_C=\log_X P_C\). Primitivity makes \(P_C\) coprime
to every \(Q_C(x)\) on its fibre. For a prime-power depth
\(q=p^j\), set
\[
 f_C(q)=\frac{\bigl|\{x\in\mathcal A_C:q\mid Q_C(x)\}\bigr|}{m_C},
 \qquad w_q=\log_X p .
\]
Include all depths of a prime appearing in the high tail when
they divide its full prime power. The depth sum of an individual
tail is then exactly its numerical weight. Put
\[
    H_s(C)=\gamma_C+\sum_q w_q\bigl(1-(1-f_C(q))^s\bigr).
\]
This is the exact expected numerical LCM weight of the target
tail and the \(s\) opposite tails chosen by independent uniform
sampling on the original fibre.

**PL7.** Fix \(\zeta>0\). The number of complete targets
with \(m_C\ge n^u\) and \(H_s(C)\ge 1-\theta+\zeta\)
is at most
\[
          O_{C_0,r,s,R,u,\zeta,L,\varepsilon}
             \bigl(n^{s(a-u)+r-\sigma+2/R}X^\varepsilon\bigr).
                                                               \tag{PL7}
\]
In the expectation defining \(H_s(C)\) the sampled LCM
weight is at most \(s+1\). Thus a fixed positive proportion
of stars have weight at least \(1-\theta\). Removing
repeated neighbours costs \(O_s(1/m_C)\). For the fixed
ordered source tuple, use (PL5) on the affine values
\(W,W-x_1,\ldots,W-x_s\), with \(W=C\). Summing over tuples
and dividing by \(n^{su}\) yields the stated exponent.

The two moments in (PL6) satisfy
\(\omega_C=\sum_q w_q f_C(q)\) and
\(\xi_C=\sum_q w_q f_C(q)^2\). When \(\omega_C>0\),
the exact lower bound in terms of these two moments is
\[
 H_s(C)\ \ge\
 \gamma_C+\frac{\omega_C^2}{\xi_C}
       \left[1-\left(1-\frac{\xi_C}{\omega_C}\right)^s\right].
                                                               \tag{PL8}
\]
Indeed the function
\(g_s(t)=(1-(1-t)^s)/t=\sum_{j=0}^{s-1}(1-t)^j\),
extended continuously at zero, is convex on \([0,1]\).
Apply Jensen's inequality with weights \(w_q f_C(q)/\omega_C\).
Equality holds when all positive depth frequencies are equal.

For example assume \(\gamma_C\ge\gamma\),
\(H=1-\theta-\gamma>0\), and \(\omega_C\ge\omega>0\).
Fix \(0<\delta<1\). If
\[
                   \xi_C\le (1-\delta)\frac{\omega_C^2}{H},
                                                               \tag{PL9}
\]
then (PL7) applies with any integer \(s\ge 2H/(\delta\omega)\).
To obtain an explicit margin, the inequality
\(1-(1-t)^s\ge st/(1+(s-1)t)\) on \([0,1]\) and
weighted Cauchy–Schwarz give
\[
 H_s(C)-\gamma_C
  \ \ge\ \frac{s\omega_C^2}{\omega_C+(s-1)\xi_C}
  \ \ge\ \frac{H}{1-\delta/2}.
\]
The last bound uses
\(1/(s\omega_C)+(1-1/s)\xi_C/\omega_C^2
\le (1-\delta/2)/H\).
Thus one may take
\(\zeta=H\delta/(2-\delta)\).
When \(H\le0\), every star qualifies from the target
tail alone, giving the same exponent by the deterministic
version of the double count.

If \(c>s(a-u)+r-\sigma+2/R\) with a fixed margin, the
cohort in (PL9) has a power-smaller count than \(n^c\).
This compares the mean common opposite depth to the
residual room left by the actual target's high tail. In
the near-complete dictionary limit \(a-u\to0\), one can
choose any fixed \(\delta>0\) before the height and
budget limits, holding \(c,\gamma,\omega\) fixed and
choosing \(r\) small and \(R\) large.

The limiting threshold \(\xi=\omega^2/H\) is sharp as
an implication from those two moments alone. In an
abstract occupancy model, use \(M\) independent depth
blocks, each of numerical weight \(H/M\) and frequency
\(\omega/H\) with \(0<\omega<H\). Then
\(\xi=\omega^2/H\) and the union never weighs more than
\(H\). For fixed \(s\) the probability of reaching
exactly \(H\) is
\([1-(1-\omega/H)^s]^M\), which tends to zero as
\(M\to\infty\). This does not assert the existence of
arithmetic fibres with those block statistics.
