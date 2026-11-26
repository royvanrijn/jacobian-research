# Nonprincipal signed modes for modular target popularity

At the prime scale used to test popular target radicals, finitely many
exceptional multiplicative modes give fewer targets than the critical
counting scale. The residual estimate concerns a signed correlation
on the actual target dictionary. Completion over all finite-field
targets cannot give the needed second moment, even on factorized
boxes obeying the three radical caps. A real near-equality cutoff
allows the exceptional modes to be removed at a smaller prime scale
when each summand dictionary is sparse in short physical intervals.

## Fixed target expansion

Use three integer dictionaries \(A,B,C\) with values at most \(N^D\),
where \(D\) is fixed, and individual radicals at most \(N\). Take
\(|A|,|B|,|C|\le N^{1+o(1)}\). Summand labels may have nonnegative
weights of uniformly subpower size. Combine repeated presentations of
an integer into one such weight. Target labels are taken distinct as
integers and counted without weight. Write \(D_A,D_B\le N^{1+o(1)}\)
for the total summand weights. In this note primes \(q\) are
sampled uniformly from \([Q,2Q]\), where \(Q=N^\beta\) and
\(\beta>1\) is fixed. Every dictionary value is therefore a unit
modulo \(q\).

Write \(I_q(c)\) for the weighted pair count with \(a+b\equiv c
\pmod q\). For unit residues \(\hat A_q(\chi)=\sum_a
w_A(a)\overline{\chi(a)}\), and similarly for \(B\). Extend each
multiplicative character by zero at zero, including the principal
one. The fixed-target Jacobi expansion is
\[
 I_q(c)=\frac{1}{(q-1)^2}
  \sum_{\chi,\psi}J_q(\chi,\psi)
      \hat A_q(\chi)\hat B_q(\psi)(\chi\psi)(c).
                                                               \tag{PM1}
\]
For each letter exclude a prescribed set of at most \(K\) characters,
including the principal character, with \(K\) fixed. The sets can
depend on \(q\). Let \(G_q(c)\) be the sum in (PM1) whose three
character coordinates avoid their respective excluded sets. This
is a signed sum. The use of \(\operatorname{Re}G_q\) below permits
excluded sets that are not closed under conjugation. In an anatomic
box one may exclude all characters annihilated by any original
fixed power on the corresponding letter.
[Jacobi expansions and bounded-order modes](ABC_MODULAR_JACOBI_INCIDENCE.md#jacobi-sums-and-the-removal-of-bounded-order-modes)

There are two exact resummations for the omitted terms. For a
fixed \(\chi\) on \(A\), the sum over all \(\psi\) is
\[
      \frac{\hat A_q(\chi)}{q-1}
                  \sum_b w_B(b)\chi(c-b).                      \tag{PM2}
\]
For a fixed target character \(\rho=\chi\psi\), summing over
\(\chi\) gives
\[
       \frac{\rho(c)}{q-1}\sum_{a,b}
                           w_A(a)w_B(b)\overline{\rho(a+b)}.    \tag{PM3}
\]
The first identity follows from multiplicative orthogonality in
the Jacobi sum. The second is the \(\rho\)-projection of the
indicator of the residue \(a+b\) on the unit group. Both hold
pair by pair, so a condition coupling \(a,b,c\) can be inserted
inside their pair contribution. When two of the three character
coordinates are specified, the third is determined. The finitely
many inclusion-exclusion intersections of omitted sets are
bounded using \(|J_q(\chi,\psi)|\le\sqrt q\), except that
the all-principal term has \(J_q(1,1)=q-2\) and is kept
separately.
[Resummation in one character](ABC_ABSOLUTE_JACOBI_SHARPNESS.md#exact-resummation-in-one-character-coordinate)

## Removal at the popularity scale

Suppose \(0\le u\le1/2\) and \(\beta>2-u\). Then there exist
\(\eta,\sigma>0\), depending only on the fixed exponent margins,
such that
\[
 \mathbb E_q\#\{c\in C: |I_q(c)-G_q(c)|>N^{u-\eta}\}
                                 \ll N^{1-u-\sigma}.           \tag{PM4}
\]
In particular one may bound the number of targets with
\(I_q(c)\ge N^u\) by those with
\(\operatorname{Re}G_q(c)\ge N^u/2\) and an exceptional set of
the mean size in (PM4). One can lower the starting threshold
by a sufficiently small fixed exponent allowance.

For the proof, write \(E_B(q)=\sum_{x\bmod q}w_B(x)^2\),
where weights at congruent labels have been added. For a
nonprincipal \(\chi\),
\[
 \sum_{c\in\mathbb F_q}\left|\sum_b w_B(b)\chi(c-b)\right|^2
                                  =qE_B(q)-D_B^2.              \tag{PM5}
\]
In fact the character correlation of two distinct shifts is
\(-1\), by the fractional-linear substitution
\((c-b)/(c-b')\), while the diagonal is \(q-1\).
For complex weights the same identity holds with \(|D_B|^2\)
and squared absolute masses. Each pair of distinct integer
labels can collide at only \(O_{D,\beta}(1)\) of the test
primes. The prime number theorem gives
\[
       \mathbb E_q E_B(q)\le N^{1+o(1)}+N^{2-\beta+o(1)}
\]
and the same statement for \(A\). On average
\(O(N^{2-\beta+o(1)})\) entries of \(C\) need to be removed
to make its reduction injective. This is already below
\(N^{1-u}\) since \(\beta>2-u\ge1+u\).

The total contribution of nonprincipal omitted modes on \(A\)
in (PM2) has complete squared norm at most
\(O_K(D_A^2E_B(q)/q)\). Average (PM5), and use the injective
target reduction and Chebyshev. At threshold \(N^{u-\eta}\)
this bounds the mean number of bad targets by
\(N^{3-\beta-2u+2\eta+o(1)}\). The same argument applies to
\(B\). Their principal modes and the all-principal term are
at most \(O(D_A D_B/q)\) pointwise.

Discard the primes for which either \(E_A(q)\) or \(E_B(q)\)
exceeds \(N^{1+\eta}\). The integer diagonal energies are
\(N^{1+o(1)}\), so the off-diagonal collision bound makes the
mean cost in targets \(O(N^{2-\beta-\eta+o(1)})\).
On the remaining primes Cauchy--Schwarz and multiplicative
Parseval give the pointwise bound
\[
 \left|\frac{1}{(q-1)^2}\sum_\chi
   J_q(\chi,\rho/\chi)\hat A_q(\chi)\hat B_q(\rho/\chi)\right|
                \ll q^{-1/2}\sqrt{E_A(q)E_B(q)} .
\]
It is at most \(N^{1+\eta-\beta/2+o(1)}\), below
\(N^{u-\eta}\) when \(\eta\) is chosen sufficiently small.
The intersections of excluded sets apart from the
all-principal term have size \(O_K(N^{2-3\beta/2+o(1)})\)
by the Jacobi bound, and are
also below the threshold. Taking \(\eta\) smaller than the
strict exponent margins proves (PM4).
[Prime-averaged collision bound (MF3)](ABC_MODULAR_JACOBI_INCIDENCE.md#jacobi-sums-and-the-removal-of-bounded-order-modes)

This supplies a reduction for any fixed number of independent
single-prime blocks in (MP1). On a block where a target radical
\(R\) has \(N^u\) admitted support pairs, choose an eligible
representative triple for each pair. There are \(N^{o(1)}\)
complete target representatives having radical \(R\), uniformly
at polynomial height. At least one \(c\) therefore has
\(N^{u-o(1)}\) distinct pairs in the Cartesian integer count.
One may use an upper envelope by three dictionaries for a
more restrictive coupled support test. Apply (PM4) with a
smaller threshold exponent and a union bound over the fixed
number of blocks. Outside a set of expected size
\(O(N^{1-u-\sigma'})\), every such target radical must have
a representative with a large positive real part of the signed
high-mode sum on each block. The representative can vary by
block. A sufficient missing estimate on a Cartesian upper
envelope would be the target-restricted signed moment
\[
      \mathbb E_q\sum_{c\in C}
                 \bigl(\operatorname{Re}G_q(c)\bigr)_+^2
                       \ll N^{1+u-\zeta}
\]
for a fixed \(\zeta>0\), with room for the threshold allowance.
Here \(t_+=\max(t,0)\). No such bound is inferred here.
[Independent blocks and support thresholds (MP1)](ABC_MODULAR_TARGET_POPULARITY.md#support-multiplicities-and-independent-blocks)

## Near-equality windows with short-interval dictionaries

Now count only pairs with \(|a+b-c|\le W\) in the modular
test, for a fixed \(W\ge1\) that may depend on \(N\).
Use the resulting truncated pair contributions to define
\(I_{q,W}\) and \(G_{q,W}\) in (PM1). Suppose in addition that
each summand dictionary has total weight at most
\(N^{a_0+o(1)}\) in every physical interval of length \(W\),
for a fixed \(0\le a_0\le1\). For \(0\le u\le1/2\), the
conclusion (PM4) for \(I_{q,W}-G_{q,W}\) holds whenever
\[
                   \beta>\max(1+a_0-u,\;1+a_0/2).
                                                               \tag{PM6}
\]
An exactly fixed discrepancy interval is used in the proof.
A fixed finite union of intervals can also be treated by
summing. The target test may require a stricter discrepancy
condition by using this window as an upper envelope.

If the window contains the whole possible discrepancy range,
use its constant indicator. Otherwise choose an integer \(P\)
of polynomial height larger than twice the sum of \(W\) and
the maximal absolute discrepancy. Use the additive expansion
of the window modulo \(P\). The normalized
coefficients \(k_t\) have sum of absolute values \(O(\log P)\),
by the geometric sum bound. Their terms insert
\(e_P(ta)e_P(tb)e_P(-tc)\), where
\(e_P(x)=\exp(2\pi i x/P)\).

Partition the three dictionaries into intervals of length
\(W\), indexed by \(i,j,k\) respectively. If
\(|a+b-c|\le W\) then \(i+j=k+O(1)\), with at most an
absolute number of offsets. Keep only these bin pairs in
the expansions since the discrepancy filter vanishes on all
others. For a fixed omitted nonprincipal \(\chi\) on \(A\),
a fixed \(t\), and one offset, the contribution for target
\(c\) in bin \(k\), after factoring out its phase, is
\[
 \frac1{q-1}\sum_{i+j=k+\text{offset}}
          \hat A_{i,q,t}(\chi)\,S_{j,q,t}(\chi;c),\qquad
 S_{j,q,t}(\chi;c)=\sum_{b\in B_j}w_B(b)e_P(tb)\chi(c-b).
\]
Here \(\hat A_{i,q,t}\) uses the phase-weighted \(A\)-bin.
By Cauchy--Schwarz its absolute value has a majorant
depending only on \(c\bmod q\):
\[
       \frac{(\sum_i D_{A_i}^2)^{1/2}}{q-1}
                 \left(\sum_j |S_{j,q,t}(\chi;c)|^2\right)^{1/2}.
\]
The same formula (PM5) for complex weights bounds its squared
norm on all residues by
\(O((\max_i D_{A_i})D_A E_B(q)/q)\). Within-bin residue
energies sum to at most the energy of the full nonnegative
summand dictionary. The triangle inequality in the squared-integrable norm,
with the additive Fourier coefficients and the finitely many
offsets and characters, preserves this bound up to subpower
factors. Apply the symmetric argument on \(B\). Thus the
omitted nonprincipal summand terms have a joint residue
majorant whose set above threshold \(N^{u-\eta}\) has
mean size at most
\[
                         N^{2+a_0-\beta-2u+2\eta+o(1)}.
\]
To pull this set back to distinct integer targets, write
\(K_q\) for its residue size and \(E_C^{\rm off}(q)\) for the
unweighted count of ordered distinct target labels congruent
modulo \(q\). The pullback has at most
\(K_q+\sqrt{K_q E_C^{\rm off}(q)}\) labels, by
Cauchy--Schwarz on the excess multiplicities in its residues.
On average the two terms cost respectively
\(N^{2+a_0-\beta-2u+2\eta+o(1)}\) and
\(N^{2+a_0/2-\beta-u+\eta+o(1)}\).
The two inequalities in (PM6) make both smaller than the
critical target count.

For the principal summand terms use (PM2) on each bin
rectangle, or its pair kernel directly. The total compatible
pair mass for every \(c\) is at most \(O(N^{1+a_0+o(1)})\).
This bounds the principal terms by
\(O(N^{1+a_0-\beta+o(1)})\). The same pointwise bound follows
from (PM3) for every omitted target mode, since each pair
has kernel at most \(1/(q-1)\) in absolute value. The all-principal term has the same \(O(N^{1+a_0-\beta+o(1)})\)
bound. The other bounded intersections have bound
\(O_K(N^{1+a_0-3\beta/2+o(1)})\). These are all below
the popularity threshold by (PM6), completing the proof.

This establishes removal of the exceptional modes for the
specified real filter. It does not replace the signed
target-restricted estimate for the remaining modes.

## Unit tests at the smaller real-filtered prime scale

Consider the real filter \(|a+b-c|\le X^\theta\) from
(SI3), with \(N=X^L\), \(X/2<c\le X\), and all three
radical caps at most \(N\). Fix \(0\le u\le1/2\) and
\(\beta>1-u\). A single-prime block now uses
\(q\in[N^\beta,2N^\beta]\); the prime may be smaller than
\(N\). Admit a record only with representatives for which
\(a,b,c\) are all units modulo \(q\). Then the true-tail
comparison (SI3) still holds up to a fixed factor in the
threshold and a power-smaller expected loss of targets.
On each block the finitely many exceptional character
coordinates can be removed pointwise. The remaining condition
for popular targets is a large positive real part of the
high-character sum on the real-filtered, unit-restricted
Cartesian envelope.

For the unit restriction, fix one true triple for each
true support pair at a target \(R\). At fixed \(L,\beta\)
its product is divisible by at most a bounded number of
the test primes, so the probability it is nonunit is
\(N^{-\beta+o(1)}\). If \(m(R)\ge2H\), the probability
that fewer than \(H\) true support pairs remain is at most
\(O(N^{-\beta+o(1)})\), by counting missing pairs in
expectation and applying Markov's inequality. For finitely
many blocks the expected number of affected true targets
is \(O(N^{1-\beta+o(1)})\). This is smaller than
\(N^{1-u}\) by a fixed power because \(\beta>1-u\ge u\).
False records still obey the upper comparison in (SI3):
the unit condition only removes possible witnesses.

For each fixed complete \(c\), the number of compatible
Cartesian pairs \((a,b)\) in the real window is
\(N^{1+\eta+o(1)}\) by (SI2), for any prescribed
\(\eta>0\) if the width exponent was chosen accordingly.
Indeed fix the smaller summand and count only the larger
summand in its short interval at height comparable to \(X\).
This avoids a short-interval condition on the smaller
summand. It also bounds the pair mass with subpower weights.
For any \(q\) restrict the summand dictionaries to units
and apply (PM1) only for unit \(c\). The pair kernels in
(PM2) and (PM3) have magnitude at most \(1/(q-1)\).
Any prescribed bounded collection of omitted modes on
any one coordinate therefore contributes
\(O(N^{1+\eta-\beta+o(1)})\) at every \(c\). An intersection
costs at most this as well, since its kernel is
\(O(q^{-3/2})\) apart from the all-principal term of size
\(O(q^{-1})\). Choose \(\eta<\beta+u-1\).
All these contributions are smaller than \(N^u\) by a fixed
power, even allowing the subpower loss from selecting
one complete representative per popular target support in
a block. Thus for every admissible block above this threshold
a bound on the number of \(R\) with a large real part of
the high-character sum gives the same bound on modular
popularity. As in (MP1), the representative may change
across blocks.
[Short-interval estimate and multiplicity test (SI2, SI3)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests)

## A signed complete-residue second moment obstruction

For any fixed \(\beta>1\) there are factorized, pairwise
prime-disjoint dictionaries with the actual caps and least
target repair of the absolute Jacobi example, with no exact
integer equality, for which
\[
 \sum_{c\in\mathbb F_q^\times}|\operatorname{Re}G_q(c)|^2
       \ge N^{2-o(1)},\qquad
 \sum_{c\in\mathbb F_q^\times}
      (\operatorname{Re}G_q(c))_+^2
       \ge N^{\min(\beta,2)-o(1)} .                            \tag{PM7}
\]
Both bounds hold on a positive proportion of test primes and
concern the unfiltered \(G_q\). The sum uses all unit residues,
instead of the actual target dictionary.
[Construction of actual capped prime-product dictionaries](ABC_ABSOLUTE_JACOBI_SHARPNESS.md#sharpness-of-the-absolute-jacobi-bound-on-primitive-balanced-boxes)

Start with the disjoint prime pools of that construction, and
independently keep each prime in each summand slot with
probability one half. A summand pair is uniquely specified
by its choices of one prime in each summand slot. The target
dictionary need not be thinned. Choose the base-\(N\) root
interval constants, or rescale \(N\) by a fixed constant, to
put all three radicals at most \(N\). Before thinning, the
summand pair count is \(N^{2-o(1)}\). For most test primes,
all but \(o(1)\) of these pairs have \(a+b\not\equiv0\pmod q\),
because each positive sum has at most \(O(1)\) test prime
divisors. All three radical caps and the lack of exact sums
are retained by thinning.

Let \(K_q(c;a,b)\) be the coefficient of one summand pair in
\(G_q(c)\). By (PM2) and (PM3), the omitted terms in its
kernel are bounded by \(O_K(1/q)\) per one-coordinate mode
and \(O_K(q^{-3/2})\) per intersection other than the
all-principal coefficient of size \(O(q^{-1})\). Indeed (PM3) gives
the exact one-pair kernel for the target-coordinate mode.
Thus \(\operatorname{Re} K_q(a+b;a,b)\ge1-O_K(1/q)\) at
each nonzero sum residue. Write each Bernoulli slot variable as
\((1+\epsilon)/2\) with an independent Rademacher sign
\(\epsilon\). The component of full degree in the thinning
variables has distinct orthogonal sign monomials for each
pair. If there are \(s\) summand slots, this gives
\[
 \mathbb E_{\rm thin}\sum_c|\operatorname{Re}G_q(c)|^2
       \ge 4^{-s}\sum_{c,a,b}|\operatorname{Re}K_q(c;a,b)|^2
       \gg_{s,K} N^{2-o(1)}
\]
on the indicated primes. For each fixed \(c\), its fourth norm is at most
\(3^{s/4}\) times its second norm. To see this, write
the polynomial in one slot as \(a+\sum_i b_i\epsilon_i\),
conditioning on the other slots. Its conditional fourth
moment is at most \(3(a^2+\sum_i b_i^2)^2\) by expanding.
Apply the triangle inequality in the second norm to the
sum of squares of the coefficients, and induct over the
slots. Orthogonality gives the corresponding sum of their
second moments. A second application of the triangle inequality bounds the second norm of
\(\sum_c|\operatorname{Re}G_q(c)|^2\) by a constant times its
mean. Cauchy then gives a fixed positive probability of
attaining at least half that mean. Exponentially unlikely
thinnings that remove more than three quarters of any slot
can be discarded. Average over primes to choose one
deterministic thinning achieving the first bound on a
positive proportion of test primes.
[Decoupled thinning moment argument](ABC_ABSOLUTE_JACOBI_SHARPNESS.md#exact-resummation-in-one-character-coordinate)

To bound the positive part write \(g=\operatorname{Re}G_q\)
and \(E_\pm=\sum_c(g(c)_\pm)^2\). Its complete sum is zero
because the principal target mode was excluded. The omitted
kernels have bound \(O_K(1/q)\) per pair, so
\(g(c)\ge-O_K(D_A D_B/q)\), as \(I_q(c)\ge0\).
If \(E_+\) already supplies half of the first bound in (PM7),
it also supplies the second bound. Otherwise
\[
 E_- \ll_K \frac{D_A D_B}{q}\sum_c g(c)_-
       = O_K\!\left(\frac{D_A D_B}{q}\sum_c g(c)_+\right)
       \ll_K \frac{D_A D_B}{q}\sqrt{q E_+}.
\]
Together with \(E_-\ge N^{2-o(1)}\) this gives
\(E_+\ge N^{\beta-o(1)}\).

At \(\beta>2-u\) and \(u\le1/2\) one has
\(\min(\beta,2)>1+u\). Thus a complete-residue bound of
order \(N^{1+u-\zeta}\) fails even for the positive part
of the signed high-mode sum. The restriction to the actual
third dictionary remains necessary.
