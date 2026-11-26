# Recovering distinct least-repaired targets from primitive odds

The exact quotient of two primitive odds fixes the possible common
part of the two complete targets. Consequently a reduced ratio of
target radicals is enough to recover both rows. A Pell or Thue record
on the target roots, or a sufficiently precise real window for the
ratio of the least-repaired heights, can supply that ratio. The
recovery counts the original primitive rows with their actual
prescribed repairs.

## Reduced target radical ratios and exact odds

At height \(X\) use distinct ordered primitive rows
\[
        C_i=\sigma U_i+\tau V_i,\qquad
        (U_i,V_i)=1,\qquad U_i,V_i,C_i\in[1,X],\qquad i=1,2,
\]
with \(\sigma,\tau\in\{-1,1\}\) fixed across the pair. In particular
\(C_i\) can denote any one fixed original role, interchanging the
roles and signs finitely many ways if needed. Write
\[
  E=\frac{U_1V_2}{U_2V_1}=\frac A B,\quad (A,B)=1,\qquad
  R_i=\operatorname{rad}(C_i),\qquad \rho=\frac{R_1}{R_2}.
\]
The odds \(U_i/V_i\) determine the original row, so \(A\ne B\).
Both \(E\) and \(\rho\) are positive.

**RRP1.** For fixed \(E\) and fixed reduced \(\rho\), the original
ordered pair count is \(O_\varepsilon(X^\varepsilon)\), uniformly
in these exact data. All caps on the original letters may be
imposed during the count.

Put \(s_i=\sigma U_i\) and \(t_i=\tau V_i\), so \(s_i+t_i=C_i\).
The common complete-target divisor \(g_C=(C_1,C_2)\) divides
the cross difference \(s_1t_2-s_2t_1=s_1C_2-s_2C_1\).
Every prime of \(g_C\) is a unit on all four signed summands
by primitivity. It is therefore a unit on the factor canceled in
reducing the cross ratio. This proves
\[
                         g_C\mid |A-B|.
\]
Write \(\rho=u/v\) in lowest terms. Then \(R_1=g_Ru\),
\(R_2=g_Rv\), where \(g_R=(R_1,R_2)\) divides
\(\operatorname{rad}|A-B|\). As \(A,B\le X^2\), enumerating \(g_R\)
costs a subpower. With \(R_i\) fixed, the complete \(C_i\)'s have
a fixed support and range over \(X^{o_X(1)}\) possibilities,
uniformly in their supports. For each of these complete targets,
the identity
\[
       ((A-B)s_1-A C_1)((A-B)s_2+B C_2)=-AB C_1C_2
\]
follows directly from \(B s_1(C_2-s_2)=A s_2(C_1-s_1)\).
The right side is nonzero and has polynomial height, so there
are only divisor-many signed factorizations. Each one fixes
\(s_1,s_2\), hence both rows. This is the signed fixed-target-pair
completion from the complete odds identity.

In particular the equal-target-radical case \(\rho=1\) costs
\(X^\varepsilon\) per exact \(E\), including distinct complete
targets on the same support. It still needs a pair relation if
each target support occurs only once.
[Complete cross-ratio factorization (TR1)](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio)
· [Prescribed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Repaired target records and their price

Fix a small positive \(L\), put \(n=X^L\), and let the inner height
limit be at this \(L\). Fix \(K>0\) and choose
\(M=\lceil K/L\rceil\), taking \(L\) small enough that \(M\ge3\).
Use the genuine least repairs of the chosen complete letters,
\[
 H_i=H_M(C_i)=R_i^M W_i,\qquad
 W_i=\prod_{p\mid C_i}p^{(v_p(C_i)-M)_+}.
\]
When the individual target radical cap is \(R_i\le n^{a_i}\)
for fixed \(a_i>0\), the heights \(H_i\) are bounded by a fixed
power of \(X\), uniformly at small \(L\). In fact
\(H_i\le X\, n^{a_iM}\). The tail \(W_i\) is at most \(X\), and
each prime on its support has original target valuation at least
\(M+1\). Thus
\[
                  \operatorname{rad}(W_1W_2)
                          \le X^{2/(M+1)}
                          \le n^{2/K}.
\]
Choosing the joint squarefree support and then both exact powers
on that support gives an \(O(n^{2/K}X^\varepsilon)\) tail catalogue,
uniformly. If the tails are already prescribed, their price is
zero. This is an upper bound on actual original tails, with no
repair multiplier counted independently.

Assume the exact odds quotient \(E\) lies in a prescribed or already
proved dictionary of size at most \(n^{\gamma+o_X(1)}\).
All additional records below can be conditional on the fixed
entries, including the tails when they are used. In a branch where
odds have been counted by a rational congruence or real strip,
the dictionary can be its actual set of exact outputs. For example,
under the FRP5 relative root and coefficient height hypotheses
and a prescribed common deep full divisor \(D\mid (C_1,C_2)\),
the exact output count is \(O(n^\eta X^\varepsilon)\) for any fixed
allowance \(\eta>0\). This step uses only that \(E\equiv1\bmod D\),
by original primitivity, and works with different complete
targets. If the chosen target role is numerically short in both
rows, with a fixed gap to their two long free roles,
the corresponding real-strip count applies. Any actual divisor
dictionary is paid before using such an output count.
[Distinct rational outputs with unit congruences](ABC_SIGNED_VALUATION_PHASES.md#a-rational-congruence-bound-for-distinct-exact-outputs)
· [Real and modular relative odds counts](ABC_FIXED_FIBRE_RELATIVE_ODDS_PAIRS.md#moving-targets-with-a-common-full-divisor)

### A reduced target norm or a nonlinear repaired difference factor

**RRP2.** After fixing \(E,W_1,W_2\), let \(G\) be an admissible
integer binary homogeneous form from a prescribed finite list.
Its coefficients have polynomial height and its degree is bounded
at each fixed \(L\). Require \(G\) to be a nondegenerate quadratic
or have at least three distinct complex projective linear factors.
Let the nonzero exact value \(G(R_1,R_2)\) have radical at most
\(n^\beta\) outside a prescribed prime pool of polynomial-height
product, where the pool may depend on all the already fixed data.
Then the pair count is
\[
                       O(n^{\gamma+2/K+\beta}X^\varepsilon).
\]
A prescribed list of forms of size \(n^{c+o_X(1)}\) may be
allowed by adding \(c\) to the exponent. If the chosen form does
not depend on the tails, the same argument requires no tail
records and removes \(2/K\). The same bounds hold with the form
evaluated on the reduced target radical pair \((u,v)\).

Indeed enumerate the outside radical of \(T=G(R_1,R_2)\),
the exact prime powers on its union with the fixed pool, and
the sign. Since \(|T|\) has polynomial height at fixed \(L\),
the cost is \(O(n^\beta X^\varepsilon)\). For unreduced inputs choose
the gcd \(g_R\) at divisor cost from \(|A-B|\). Writing \(d=\deg G\),
the remaining coprime pair satisfies \(G(u,v)=T/g_R^d\).
For a form with at least three distinct projective factors,
Evertse's coefficient-uniform bound for primitive solutions of
the Thue--Mahler equation, using the primes of the nonzero
right side, gives height-subpower completions. For a
nondegenerate quadratic \(G(u,v)=a u^2+buv+c v^2\) with \(a\ne0\),
complete the square as
\[
   (2au+bv)^2+(4ac-b^2)v^2=4aT/g_R^2 .
\]
Both coefficients on the left and the right side are nonzero.
All heights are polynomial, so the Vaughan--Wooley bound gives
height-subpower completions, including negative or square
discriminants. With \(a=0\), the split equation \(v(bu+cv)=T/g_R^2\)
gives a divisor count. RRP1 finishes the recovery for each
reduced ratio. In the reduced-input version apply the same
arguments directly at level \(T\). No common multiplier is
summed over beyond the divisor bound and the original
supported powers.
[Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14)
· [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)

An example is any primitive irreducible nonlinear factor of
\(W_1U^M-W_2V^M\), evaluated at \(R_1,R_2\), when \(H_1\ne H_2\).
The binomial has distinct projective roots. Since \(M\ge3\),
it has a nonreal root and thus a nonlinear factor over
\(\mathbb Q\). Every such quadratic is nonsingular and higher
degree factors have at least three distinct projective roots.
The factor coefficients have polynomial height at fixed \(M\).
One can also take the entire binary form, to count by the
complete repaired difference. Equal repaired heights have
equal radicals and are covered by RRP1 without the tails.
Unlike the normalized repair-lift recovery, the fixed odds
already determines the relation between the original free
letters. The only binary completion is on the target radicals.
In particular there is no separate charge for a reduced
summand ratio or a second power-defect record.
[Repair lifting with summand and target records (WD3)](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#pair-records-and-their-count)

### Short real windows for least-repaired target ratios

**RRP3.** After fixing \(E,W_1,W_2\), suppose \(H_1/H_2\)
is in a relative window
\[
                  \left|(H_1/H_2)/C_*-1\right|\le \delta
\]
with \(0\le\delta\le1/2\) about a prescribed positive real center
\(C_*\). This record has pair count
\[
                  O_\varepsilon\bigl((1+\delta n^{a_1+a_2}/M)
                                      X^\varepsilon\bigr).
\]
The bound is uniform in \(C_*\). If the conditional center
catalogue has size at most \(n^{b+o_X(1)}\) and
\(\delta\le n^{-t}\), the full pair exponent is at most
\[
                   \gamma+2/K+b+(a_1+a_2-t)_+
\]
up to an arbitrary small height allowance in normalized
units. A fixed power precision \(\delta\ll X^{-\lambda_*}\)
with \(\lambda_*>0\) gives exponent \(\gamma+2/K+b\)
for all sufficiently small fixed \(L\).

Since \(H_i=R_i^M W_i\), the window puts
\(\log(R_1/R_2)\) in an interval of length \(O(\delta/M)\).
If two distinct reduced ratios \(u/v,u'/v'\) in the
allowed box have \(u/v>u'/v'\), their quotient differs
from one by
\((uv'-u'v)/(u'v)\ge n^{-a_1-a_2}\).
Their logarithmic separation is at least a fixed
multiple of \(n^{-a_1-a_2}\), uniformly in the numerical
position of the interval. This gives the claimed ratio
count, then RRP1 gives pair recovery. If every
eligible pair has \(g_R\ge n^g\), first enumerate the
possible \(g_R\)'s from \(|A-B|\) at divisor cost. The
height caps on the reduced inputs then give the
stronger \(a_1+a_2-2g\) in place of \(a_1+a_2\).
For example any common divisor \(F\) of the original
targets yields \(g\ge\log_n\operatorname{rad}(F)\).
The record is a real window on repaired heights
\(H_i\) about a prescribed center. The proof allows
unequal complete targets with their fixed actual
least-repair tails.
This is a different hypothesis from the residue
of a corrected repair ratio on a two-summand lock.
[Summand-lock residue version (TL3)](ABC_TWO_TARGET_SUMMAND_LOCKS.md#repair-multiplier-residues-and-target-support-recovery)

In either arm, if the exponent lies below \(1/3\) by a
fixed margin, the primitive rows incident to the
eligible pairs also have a saving, as each pair has
only two endpoints. The pair estimates thus apply
with distinct original targets and a single
representation per target support. An unrestricted
application still needs a relation reaching enough
actual three-cap rows with the stated conditional
records. A real window on the unrepaired \(C_i\)'s
alone does not imply the window on \(H_i\)'s, and the
least-repair operation does not by itself bound the
radical of a difference factor.


## A necessary large radical on normalized repaired differences

The preceding recovery charges a cheap target-difference
record. There is a general obstruction to extracting
such a record from typical pairs by choosing a repair
order with many fixed divisors. This obstruction needs
no compression of primitive odds.

**RRP4.** Fix \(d>u\ge0\) and let an original family have at
least \(n^d\) rows and at most \(n^{u+o_X(1)}\) rows at
each individual target radical \(R=\operatorname{rad}(C)\).
For any prescribed \(B>0\) one can choose a fixed integer
\(Q\ge2\), a fixed \(K>0\), and a positive \(\sigma\),
with the following property. At sufficiently small fixed
\(L\) use
\[
                       M=Q\left\lceil\frac{K}{QL}\right\rceil .
\]
Given each reference row a prescribed polynomial-height
prime pool \(P\) (possibly different for every reference),
all but a proportion \(O(n^{-\sigma})\) of the other
rows have distinct
radical \(R'\), and the outside-\(P\) radical of
\[
                  \frac{H_M(C')/W'-H_M(C)/W}{(R,R')^M}
                  = (R'/(R,R'))^M-(R/(R,R'))^M
\]
exceeds \(n^B\). One can also pass to a subfamily with a
single actual tail \(W_0\) at a cost of at most
\(n^{1/K+o_X(1)}\) in cardinality. After increasing \(K\)
in terms of \(d-u\), in that subfamily the same conclusion
holds for the normalized difference of the actual
repaired heights:
\[
                 \frac{H_M(C')-H_M(C)}
                             {\gcd(H_M(C'),H_M(C))}.
\]
The target role, the height cutoff and all
three original radical restrictions remain those of the
original family. The obstruction is on full-pair
extraction of low radical differences and does not
exclude all rare low-cost incidences.

To prove this, choose a real \(v\in(0,d-u)\), for example
\(v=(d-u)/4\), and fix \(Q\) with \(\tau(Q)v>B\). Choose
\(K\) large enough that \(d-1/K-u>v\) by a fixed margin.
At a fixed budget \(M\) is finite, \(Q\mid M\), and
\(ML\) lies in \([K,K+QL]\). The stated tail count costs
at most \(n^{1/K+o_X(1)}\) for a single \(W\). Choose the
largest tail class if using the normalized actual
difference; its size is at least \(n^{d-1/K-o_X(1)}\).
Without that choice the original row lower bound is
\(n^d\).

Now fix \(R\) from the reference, and choose
\(g=(R,R')\mid R\) at height-subpower cost. Put
\(r=R/g,\ s=R'/g\), so \((r,s)=1\). For each fixed
\(e\mid Q\), consider the homogeneous cyclotomic
value \(\Phi_e(s,r)\). When \(R'\ne R\) it is
nonzero. Enlarge the prescribed pool by the primes
of \(Q\). The signed values of outside radical at
most \(n^v\) have at most \(n^{v+o_X(1)}\) possible
values, uniformly in the prescribed pool. Their
height is polynomial at fixed \(Q\). With \(r\)
fixed, a nonzero value of the monic polynomial
\(\Phi_e(S,r)\) allows at most \(\varphi(e)\) possible
\(s\). Each resulting radical accounts for at most
\(n^{u+o_X(1)}\) original rows. Hence all the exceptional
neighbors for the finite list of divisors have at
most \(n^{u+v+o_X(1)}\) rows, and the coincident
target radical adds at most \(n^{u+o_X(1)}\).

For the remaining neighbors, outside the enlarged
pool the \(\Phi_e(s,r)\) have disjoint support.
Indeed at primes not dividing \(Q\) the projective
zeros of the factors of \(S^Q-r^Q\) do not coincide,
by \((r,s)=1\) and separability modulo the prime.
As \(Q\mid M\), the product of the selected factors
divides \(s^M-r^M\). Their combined outside radical
exceeds \(n^{\tau(Q)v}>n^B\); this also holds outside
the initial smaller pool. In the common-tail class
\(\gcd(H_M(C'),H_M(C))=W_0 g^M\), so the normalized
actual difference is the same integer. The strict
margins above give the uniform power-fraction bound
at sufficiently large inner heights.

For a balanced near-full family at original exponent
\(1/3-o(1)\), the fixed-letter count and subpower
enumeration of complete values on each support
give \(u\le1/6+o(1)\). Thus the gap needed by RRP4
is available throughout the balanced target
multiplicity range, including occurrence. This
estimate concerns the normalized least-repaired
difference and is separate from the physical
difference estimates on the original complete
values.
[Original fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Physical shift radical estimates (FE1–FE3)](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md)

The factor argument in fact gives a bound for any
integer \(M\) fixed at the inner height limit. For
every fixed positive \(v<d-u\), at most
\(n^{u+v+o_X(1)}\) neighbors have outside radical at
most \(n^v\) on any one of the finitely many
\(\Phi_e(s,r)\), \(e\mid M\), after excluding the
primes of \(M\). The outside radical of the whole
normalized power difference is consequently greater
than \(n^{v\tau(M)}\) for the other neighbors. For
\(M\) varying only with the outer budget the fixed
inner constants may depend on its actual value.
The fixed \(Q\) formulation above keeps the chosen
factor list and the lower-bound margin fixed as
the budget decreases.
