# Reversed product signatures on fixed and moving primitive targets

Fix \(k\ge2\), \(n=X^L\), and a family \(\mathcal A\) of positive primitive triples \(a_i+b_i=c_i\), \(a_i<b_i\) and \(X/2<c_i\le X\). All individual raw radical caps are bounded in base \(n\). Put
\[
 P=\prod_{i<k}a_i,\quad Q=\prod_{i<k}b_i,\quad (a,b)=(a_k,b_k),\quad
 t=(-1)^{k+1}.
\]
Use the two integer signatures
\[
 S=Qb+tPa,\qquad S'=Qa+tPb,
\]
with \(g=\gcd(Qb,Pa)\), \(g'=\gcd(Qa,Pb)\), \(D=S/g\), \(D'=S'/g'\). The signed second signature is recorded in its exact signed value. The first is positive since \(P<Q\) and \(a<b\). If \(S'=0\), necessarily \(t=-1\), and the other \(k-1\) rows determine the last odds \(a/b=P/Q\), hence determine the last primitive row. There are at most \(|\mathcal A|^{k-1}\) such ordered tuples, or at most \(N_c^{k-1}\) among tuples with all targets equal to \(c\). In what follows both signatures are nonzero.

## PSN1. Completion with no anchor row

For all choices of original targets, the records \(D,D',g,g'\) determine the complete ordered input tuple at \(O_{k,\varepsilon}(X^\varepsilon)\) cost. The exact unreduced integers give
\[
 S+S'=c_k(Q+tP),\qquad S-S'=(b-a)(Q-tP).
\]
All four factors on these right sides are strictly positive. In particular neither product vanishes. Choose their ordered divisors to obtain \(c_k,Q+tP,b-a,Q-tP\). This recovers \(a,b,P,Q\) exactly. Factoring \(P\) into its \(k-1\) labeled positive summands \(a_i\) and \(Q\) into the \(k-1\) labeled \(b_i\)'s costs a height-subpower number of choices. Each \(c_i=a_i+b_i\) is now determined, and all original conditions can be imposed. All integers in the proof have height \(X^{O_k(1)}\). Fixing \(|DD'|\) gives \(D,D'\) up to divisor choices and the sign of the second. Thus one exact product records the two reduced signatures at subpower multiplicity.

At a fixed equal original target \(c\), every input summand is a unit modulo \(c\), and each \(a_i\equiv-b_i\bmod c\). The choice of \(t\) makes \(S,S'\equiv0\bmod c\); both gcds are units there. Hence \(c\mid D,D'\) in the equal-target sector. For any prescribed prime pool \(P_0\le X^B\) and numbers \(\beta,\mu\ge0\), suppose on one equal-target fibre the joint radical of \(|DD'|\) outside \(cP_0\) has cap \(n^\beta\) and the actual gcd record \((g,g')\) has at most \(n^\mu X^\varepsilon\) possibilities. Enumerate exact \(D,D'\) by choosing their outside squarefree product and then all powers on the resulting polynomial-height support. The completion gives the count
\[
 O(n^{\beta+\mu}X^\varepsilon).
 \tag{A}
\]
This is uniform in \(c\) and the prescribed pool. If the equal-target fibre has \(N_c\ge n^{u-\theta}\) then any fixed gap \(\beta+\mu<ku-k\theta\) gives a power-small proportion of its tuple population, with the exponent allowances small enough. The private-role application of FL1 below supplies a generic \(\mu\le k(r_a+r_b-4u+2\xi)/2\), because a prime of either gcd occurs in at least two different original slots, whether or not the last slot was reversed. Both gcds may be counted jointly for this bound. The asymmetric one-sided improvement for an unreversed rooted system is not implied here. In the common-target limit \(u\to1/6\), \(r_a,r_b\to1/3\), the joint outside fresh mass is therefore at least \(k/6-o(1)\). At a fixed repaired target the original-target partition and Hölder argument below retain the full tuple exponent, and (A) is summed over the height-subpower many original complete targets.

## PSN2. Outside support for moving targets

For an unrelated family of moving targets, the same unconditional completion gives \(O(n^{\beta+\mu}X^\varepsilon)\) tuples with joint output radical outside a *prescribed* pool \(P_0\le X^B\) capped by \(n^\beta\), and actual gcd catalogue price \(\mu\). There is no complete target dictionary cost in this assertion. It counts only genuine original rows, and hence only their primitive quotients when grouped by least repair.

Consider a height-first violating sequence for the balanced leading coefficient, with budget tending to zero after arbitrarily large heights at each fixed budget, family exponent approaching \(1/3\), and all three original caps \(r_a,r_b,r_c\to1/3\) in base \(n\). Let \(\xi>0\) be fixed. On all but a power-small proportion of independent ordered \(k\)-tuples, the joint gcd radical has normalized mass at most \(O_k(\xi)\). Condition on every row but slot \(i\); its overlap with the union of supports of the other rows is a prescribed polynomial-height pool. EB5 gives a power saving if at least \(\xi\) of original raw mass is on that pool, since the global original total budget is \(1+o(1)\). If the original caps are stated in base \(X\) and only give \(\operatorname{rad}(abc)\le n^{1+\vartheta}\), apply EB5 at fixed dilated budget \((1+\vartheta+\vartheta')L\) with fixed positive \(\vartheta'\) absorbing the factor-two shell variation. Take \(\vartheta,\vartheta'\) sufficiently small compared with \(\xi\) before the small-budget threshold. A pool supported on one prescribed polynomial-height integer permits prime price zero in EB5, and each captured \(\xi\log n\) lowers the normalized count exponent by \(\xi/3\) before the allowances. A prime of either gcd contributes to at least two such overlaps. A fixed \(k\) union bound gives the claimed generic estimate. Therefore the joint output radical outside any prescribed pool is at least \(k/3-o(1)\) in the limiting sense on a generic tuple population.

One can additionally remove the varying original target supports \(R=\operatorname{rad}(\prod_i c_i)\) at negligible normalized raw cost. For each slot \(j\), each unreduced output \(S\) or \(S'\) is \(L_j f_j+M_j h_j\), where \(\{f_j,h_j\}=\{a_j,b_j\}\) and \(L_j,M_j\) are signed complementary monomials in the other rows. For \(p\mid c_j\) its reduction is \((M_j-L_j)h_j\), with \(h_j\) a unit. Thus a prime of \(c_j\) on the output is in the fixed polynomial-height integer \(M_j-L_j\), unless that coefficient difference is zero. For any vanishing test, choose another slot \(i\ne j\) and fix all rows except \(i\): the equality of complementary monomials (one involving \(a_i\) and the other \(b_i\)) fixes the exact positive odds of that row or is impossible. There are \(O_k(|\mathcal A|^{k-1})\) tuples failing any such test. For a nonzero test, condition on the other rows and apply EB5 to row \(j\) with this fixed pool. At the balanced endpoint all but a power-small proportion of tuples have at most \(O_k(\xi)\) total raw mass on the intersection of the outputs and the varying original target supports. In particular restricting \(\operatorname{rad}(|DD'|)\) outside \(R P_0\) to size \(n^{k/3-\Delta}\) for a fixed \(\Delta>0\) excludes all but a power-small proportion at sufficiently late budgets. Indeed on good tuples their radical outside the fixed \(P_0\) is at most \(n^{k/3-\Delta+O_k(\xi)}\), so apply the completion and take \(\xi\) sufficiently small. This gives a necessary arithmetic constraint on full-size tuple populations.

## PSN3. The actual gcd price in a fixed target fibre

Let \(u>0\) and \(\xi>0\) be fixed with \(\xi<2u\), and let \(0<\theta<\xi/8\). Suppose \(N_c\ge n^{u-\theta}\). For sufficiently small fixed \(L\), and all sufficiently large \(X\) at this fixed \(L\), all but \(O(N_c^k n^{-\xi/8})\) tuples satisfy, for every slot and each of its two free roles,
\[
 \sum_{\substack{p\mid w_i\ p\nmid P_0\prod_{j\ne i}a_jb_j}}\log_n p
                        \ge 2u-\xi,\qquad w_i=a_i\ \text{or}\ b_i.
 \tag{ON2}
\]
There is no assumption of independence within a row and no restriction on the valuation orders at the private or shared primes. Condition on all other \(k-1\) rows. Their free supports adjoined to \(P_0\) give a prescribed pool of polynomial height depending only on \(B,k\). FL1 with the fixed original complete \(c\) bounds candidates with either displayed outside cap at most \(2u-\xi\) by \(O(n^{u-\xi/2+\eta}X^\varepsilon)\). With, for example, \(\eta<\xi/16\), \(\varepsilon/L<\xi/16\), the ratio of this count to \(N_c\) is \(O(n^{-\xi/4})\). Sum over the \(2k\) slot-role pairs and weaken the margin if desired. These choices take \(\eta\) and the required \(L_0\) before the fixed-budget height limit, and take the allowed \(\varepsilon\) at each fixed \(L\). This is the same prescribed-pool private-role argument as the generic reciprocal denominator result in ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md.

A prime of any \(g_e\) occurs in at least two distinct input slots, one in each original free role, by primitivity of each row. Each of those occurrences is outside the private contribution in (ON2). Thus on good tuples
\[
 \log_n\operatorname{rad}\Bigl(\prod_e g_e\Bigr)
 \le \tfrac k2(r_a+r_b-4u+2\xi).
 \tag{ON3}
\]
Namely the total nonprivate mass across all \(2k\) free slots is at most \(k(r_a+r_b-4u+2\xi)\), and the joint gcd radical is counted twice in this total. This does not charge repeat appearances of a prime across different \(g_e\). If the right side is negative then there are no good tuples, already obsolete under the assumed lower bound on \(N_c\) and FL1. One can improve (ON3) using \(\min\big(k(r_a-2u+\xi),k(r_b-2u+\xi)\big)\), because each gcd prime is represented in at least one nonprivate first and one nonprivate second summand. Equation (ON1) also works with any better available actual conditional gcd catalogue, with no change to the proof.

For the paired outputs of PSN1, a prime of either
gcd still occurs in at least two distinct slots.
Thus the symmetric bound
\(k(r_a+r_b-4u+2\xi)/2\) applies to the union of
these two gcd radicals. Combine it with (A).
At the balanced \(u=1/6\) endpoint, the joint
fresh radical of the paired outputs is at least
\(k/6-o(1)\) outside the fixed target and the
prescribed pool. The one-sided refinement above
applies to the unreversed gcd.

There is also an unpaired cubic statement. Fix
\(c\), choose three sorted rows, and put
\(D=(a_1a_2a_3+b_1b_2b_3)/g\) with the actual
gcd \(g\). Prescribe \(D,g\) and row one,
written \(A+B=c\), \(A<B\). Writing \(x,y\) for
the first summands of rows two and three gives
\[
 ((A+B)x-Bc)((A+B)y-Bc)
                   =(A+B)gD-ABc^2.
\]
Both factors are strictly negative since \(x,y<c/2\)
and \(Bc/(A+B)>c/2\). Divisor-many choices recover
the two unknown rows. When \(N_c\ge n^{u-\theta}\),
the count with a reduced numerator cap \(n^\beta\)
outside \(cP_0\) and gcd record price \(\mu\) is
\(O(N_c n^{\beta+\mu}X^\varepsilon)\). Thus a full
balanced cubic population at \(u=1/6\) needs
fresh mass at least \(1/3-o(1)\). If the first two
rows are allowed to be opposite orientations of
one row, the reduced sum is exactly \(c\) for
every third row. The common sorting avoids this
exact exceptional configuration.

A fixed least-repaired target determines the original
target radical. There are only height-subpower many
original complete targets \(c\) on that support at
each fixed \(L,M\). In a repaired fibre of size \(N\),
if it splits into \(J=X^{o_X(1)}\) exact \(c\)-fibres,
then the number of ordered tuples wholly on one exact
\(c\) is at least \(N^k/J^{k-1}\) by Hölder.
Fibres with size below \(n^{u-\theta}\) contribute
a power-small part when \(N\ge n^{u-o(1)}\) and
\(\theta>0\) is fixed. Consequently the fixed-target
bounds apply to tuples sharing one original complete
target inside a genuine repaired fibre. These counts
use the original primitive quotients and the actual
least-repair cutoff.
[Fixed-letter radical-cap bound FL1](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Original effective support price EB5](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
