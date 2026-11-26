# Full-divisor completion for complementary-word collisions

Write each genuine positive primitive row as \(a+b=c\) and put \(x=a/c\); its complementary coordinate is \(1-x=b/c\). Fix an integer \(m\ge2\). The two exact word values for two lists \(x_1,\ldots,x_m\) and \(y_1,\ldots,y_m\) are
\[
 U=\frac{\prod_i x_i}{\prod_j y_j},\qquad
 V=\frac{\prod_i(1-x_i)}{\prod_j(1-y_j)}.
\]
All input complete values are at most \(X\). On a nonempty exact fibre the heights of \(U,V\) are at most \(X^{2m}\). The collision is nontrivial when the two multisets of original primitive rows differ.

## FC1. The complete common target and its radical

**FC1.** Suppose the positive-side rows all have the same prescribed complete target \(c\). The negative-side rows may come from any common dictionary of \(T\) genuine primitive rows, including different complete targets. For fixed exact \(U,V\), the number of nontrivial ordered pairs of lists is
\[
 O_{m,\varepsilon}(T^{m-1}X^\varepsilon)
\]
uniformly in the specified values and target. In the common-target case \(T=D_c\) can be the size of the one target fibre. This statement does not require a radical cap.

First cancel rows occurring on both sides, with as many cancellations as their smaller multiplicity. There are boundedly many patterns at fixed \(m\). The shared canceled slots may be chosen in at most \(T^s\) ways if \(s\) pairs are canceled. The remaining lists have equal length \(q=m-s\ge1\), disjoint multisets, and the same exact word outputs. It suffices to bound disjoint shorter lists by \(O(T^{q-1}X^\varepsilon)\).

Let \((A_j',B_j',C_j')\) denote negative original rows; use \(p_U,p_V\) for the reduced positive numerators of \(U,V\). The equations after clearing denominators give
\[
 a_i\mid p_U c^q\prod_j A'_j,\qquad
 b_i\mid p_V c^q\prod_j B'_j.
\]
Since each \(a_i,b_i\) is coprime to \(c\), omit \(c^q\) on both right sides. For each fixed positive slot \(i\), allocate the exponents of each prime of \(a_i\) among the known numerator \(p_U\) and the \(q\) \(A'_j\)'s, using no more than the available exponent at any source. Allocate the analogous exponents of \(b_i\) among \(p_V\) and the \(B'_j\)'s. Let \(L_{ij}\) be the combined natural logarithmic weight assigned to negative slot \(j\) in the two allocations. Then
\[
 \sum_j L_{ij}\le\log(a_i b_i).
\]
At most one \(j\) has \(L_{ij}>\log\max(a_i,b_i)\), because \(2\log\max(a_i,b_i)\ge\log(a_i b_i)\). Among the \(q\) negative slots choose a single \(j\) that is not forbidden for any of the first \(q-1\) positive slots. Fix the \(q-1\) other negative rows (cost \(T^{q-1}\)). For each of the first \(q-1\) positive slots the same allocations show that
\[
 F_i=\gcd\!\left(a_i,p_U\prod_{k\ne j}A'_k\right),\qquad
 G_i=\gcd\!\left(b_i,p_V\prod_{k\ne j}B'_k\right)
\]
have \(F_i G_i\ge a_i b_i e^{-L_{ij}}\ge\min(a_i,b_i)\).
Their possible exact values are divisors of fixed polynomial-height records, hence have height-subpower cost. These two divisors are coprime by original primitivity. The original relation with the actual fixed \(c\) gives \(a_i\equiv0\bmod F_i\), \(a_i\equiv c\bmod G_i\). For the shorter free role, fix its half-open dyadic box \([A,2A)\). The box has length \(A\le F_iG_i\) and contains at most one member of the combined CRT class in that role. Thus the first \(q-1\) positive slots are recovered at height-subpower cost. The other positive and negative fractions \(x,y\) satisfy
\[
 x/y=u,\qquad (1-x)/(1-y)=v
\]
for known positive rationals \(u,v\). The equations are affine in \(x,y\). If \(u=v=1\), they force \(x=y\), forbidden by the disjoint lists. If \(u=v\ne1\) there is no solution, and otherwise there is at most one. The primitive row is determined by its reduced fraction. This also proves the case \(q=1\). Accounting for all cancellations and divisor choices proves FC1.

**Modular consequence.** On an FP1 chart (or an SC3 chart for the same fixed signed word) for the weights \(m\) times \(+1\) and \(m\) times \(-1\), the modular identity permits \(O(n^\eta X^\varepsilon)\) exact pairs \(U,V\), with arbitrarily small fixed allowance at sufficiently small fixed budgets. Thus the number of nontrivial equal \(m\)-fold modular products with all \(2m\) rows on a common complete target is
\[
 O\!\left(n^\eta X^\varepsilon\sum_c D_c^{m-1}\right).
\]
In every regular cohort with positive original complete-target multiplicity exponent \(u\), this gives a strict gap from the full same-target moment \(\sum_c D_c^m\). For \(m=3\) this strengthens CE2's sufficient threshold \(u>5/36\) to every fixed \(u>0\). Repetitions require no additional estimate. If counting genuine minimally repaired quotients, partition each fixed repaired target first into its height-subpower many exact original complete targets, as in the fixed-support count. The common-original-target moment loses at most that subpower factor to a fixed power by Hölder; FC1 concerns these tuples. The count is on actual primitive quotients.


**Common target radicals and different targets on the two sides.** It is enough that the positive list's complete targets all have one prescribed radical \(R_+\). Enumerate its \(q\) exact complete targets at height-subpower cost on the known support, and use them individually in the CRT congruences. In the preliminary divisibility arguments replace \(c^q\) by \(\prod_i c_i\), which is still a unit relative to every positive summand by original primitivity and equality of target supports.

Moreover, let both lists have common target radicals within their own list, and permit those radicals \(R_+,R_-\) to vary independently. If \(D_R\) denotes the full original row multiplicity of the radical \(R\) in the family, then for each exact \(U,V\) the count becomes
\[
 O\!\left(X^\varepsilon\sum_R D_R^{m-1}\right).
\]
For the disjoint shortened lengths \(q\ge2\), the \(q-1\) negative anchors determine \(R_-\). For \(p\mid R_+\), \(p\nmid R_-\), all positive numerators are units at \(p\) while all negative target valuations vanish; hence \(v_p(U)<0\). Thus \(R_+\) divides the squarefree support of \(R_-\operatorname{den}U\), and is enumerated at divisor-subpower cost after fixing the negative anchors. Proceed with the preceding proof. If canceled slots occur, both target radicals coincide; for a fixed such \(R\) they cost \(D_R^s\). The case \(q=1\) requires no anchor to complete the residual fractions by the affine equations, and the \(m-1\) canceled slots already contribute only \(\sum_R D_R^{m-1}\). This proves the bound. It includes lists with fixed least-repaired targets separately on the two sides, since a least repair preserves the complete radical of the original marked letter. All counts are of the original primitive quotients.


There is also an image-size consequence in the
positive multiplicity range. Keep a bin of target
radicals satisfying \(D_R\ge n^u\) for a fixed \(u>0\).
Take all unordered \(m\)-multisets of rows with one
common target radical and look at their modular
product records from the same FP1 probes. Distinct
multisets having equal records give at least one
ordered nontrivial collision. The number of
such collisions over all pairs of target radicals is
\(O(n^\eta X^\varepsilon\sum_R D_R^{m-1})\).
For \(\eta,\varepsilon/L\) small compared with \(u\)
this is power-smaller than \(\sum_R D_R^m\).
Consequently the number of distinct product records
is
\[
                 (1-o(1))\sum_R \binom{D_R+m-1}{m}
\]
on these charts. This holds through the actual
least-repair restriction. It is a lower bound for
the occupied image, including images shared between
different original target radicals.

For arbitrary prescribed exact outputs the \(m-1\) exponent cannot be lowered uniformly. Choose any two distinct rows \(x,y\) in one sufficiently large complete-target fibre, fix the two output ratios to \(x/y,(1-x)/(1-y)\), and append the same \(m-1\) arbitrarily chosen rows on both sides. This gives \(D_c^{m-1}\) distinct nontrivial ordered pairs (with fixed positions for \(x,y\)). The obstruction concerns arbitrary exact output fibres, and says nothing about the sharp nontrivial exponent at the single exact identity.

## FC2. Moving targets on generic positive tuples

Let \(n=X^L\), and let \(\mathcal A\) be a balanced near-saturating cohort of \(T\) original positive primitive rows with \(X/2<c\le X\), all original radical caps at most \(n^r\) with \(r\to1/3\), and \(T=n^{1/3+o(1)}\) in the height-first regime. For the EB5 deletion one can instead assume their total original budget is \(1+o(1)\) in base \(n\). Use one common limiting profile box. By the pair nonatomic bound (GE1), some original role has a positive fixed amount of nonatomic mass. For any sufficiently small fixed \(\mu>0\), choose a compact positive band in which one can take full powers from that role of total normalized raw cost in \([\mu/2,\mu]\). This is achieved for every late-box row after removing neighborhoods of finitely many large atoms, by continuity and greedy selection with each prime cost at most \(\mu/2\). Thus each row can provide that full divisor \(F\) with
\[
 \operatorname{rad} F\le n^\mu,\qquad F\ge X^\nu
\]
for some fixed \(\nu>0\) depending on the band and \(\mu\). The choice of its exact value costs \(O(n^\mu X^\varepsilon)\); no dictionary of the host complete original role is charged.

For a constant \(R>0\), record for every positive slot the three full tails where \(L v_p>R\). Their exact catalogue costs \(O(n^{3m/R}X^\varepsilon)\), since every individual role has height at most \(X\). Prescribe a fixed overlap tolerance \(\xi>0\). Call a positive tuple eligible if its rows are pairwise distinct and for each slot the raw overlap of its whole original radical with the union of those of the other positive slots is at most \(\xi\). At fixed \(m,\xi\), EB5 shows that the ineligible tuples are a power-small proportion in the balanced near-saturating cohort: condition on all but one positive slot and apply the prescribed polynomial-height pool discount to that row with an allowance below \(\xi/3\). If the original cap is only stated in the \(X\)-base, take a sufficiently small fixed dilation in the total budget to fit the shell; the dilation and the saturation errors can be made small compared with \(\xi\) before the inner height limit. Eligibility survives deleting any positive slots.

**FC2.** For a prescribed exact \(U,V\), the number of nontrivial pairs of lists with eligible positive side is
\[
 O\!\left(T^{m-1} n^{r/2+\mu+3m/R+\eta} X^\varepsilon\right).
\]
Here \(R\) may be any sufficiently large fixed value; after choosing the band, \(\mu,\nu,R\) and a sufficiently small fixed \(\eta>0\), choose \(\xi\) sufficiently small in terms of \(\nu/R\). The budgets are then fixed sufficiently small. All constants and subpower allowances are uniform in the prescribed exact word values, prescribed divisors and coefficient heights under these bounds.

As before cancel shared rows first. Prove the result for the disjoint shortened lists of length \(q\le m\) and spend \(T^{m-q}\) on canceled rows. For \(q=1\) the two affine ratios give uniqueness. For \(q\ge2\) fix \(q-1\) negative rows, leaving one negative slot, and buy only the cheap full divisor \(F\) on the remaining negative row. Also fix the three positive tails on each of the \(q\) remaining slots. Use dyadic boxes for the original positive summands \(a_i,b_i\).

For each positive slot distinguish its tail (fully known), its non-tail primes shared with other positive slots (total numerical log-height at most \(R\xi\log X\)), and its remaining private non-tail primes. For a private prime \(p\) at a positive slot \(i\) with exponent \(e\), compare the two valuation vectors of \(U,V\). Primitivity gives the three rays \((1,0),(0,1),(-1,-1)\), for \(a,b,c\) respectively; no other positive slot has a contribution at that prime. On role \(a\), the first valuation gives \(e\le v_p(p_U)+\sum_j v_p(A'_j)\). On \(b\), the same holds with \(p_V\) and \(B'_j\). On \(c\), the first valuation gives \(e\le v_p(q_U)+\sum_jv_p(C'_j)\), where \(q_U\) is the reduced denominator of \(U\). (The opposite free-role terms on each side of these inequalities have the favorable sign.) In each role, allocate those private exponents among the indicated known numerator or denominator, the fixed negative rows, the known factor \(F\) of the unanchored negative row in its actual role, and that row's remaining parts. Let \(L_i\) be the combined base-\(X\) numerical weight that must be assigned to the unanchored row's remaining parts. Private supports are disjoint between all positive slots. Their total allocations to these remaining parts therefore satisfy
\[
 \sum_i L_i\le \log_X(A'B'C'/F)\le 3-\nu
\]
for the unanchored negative row. Choose \(\omega>0\) with \(3\omega<\nu\) and \(R\xi<\omega/2\). At most two positive slots have \(L_i>1-\omega\).

For every other positive slot the allocation gives three exact divisors \(D_{i,A}\mid a_i\), \(D_{i,B}\mid b_i\), \(D_{i,C}\mid c_i\), made from its known tails and the portions of the private primes supplied by the known records. Their product is at least
\[
 a_i b_i c_i X^{-R\xi-L_i}
       \ge \tfrac12 a_i b_i X^{\omega-R\xi}.
\]
The divisors may be enumerated at height-subpower cost: each divides the known product of the corresponding positive tail, the specified reduced numerator or denominator, the \(q-1\) fixed negative letters in its role, and \(F\) if it has the same role. The three divisors are pairwise coprime by original primitivity. Two candidate original primitive rows with the same divisors have determinant \(a_i b_i'-a_i' b_i\) divisible by all three, using \(c=a+b\) for the \(C\)-divisor. The determinant has absolute value at most \(4A_i B_i\) on the fixed dyadic boxes, which is smaller than the divisor product for every sufficiently large \(X\). The determinant vanishes, so the primitive rows agree. This recovers any \(q-2\) of the good positive slots. There are only boundedly many choices of which slots to leave when \(q\) is fixed.

The remaining two positive and one negative slots are distinct. Positive eligibility gives distinctness of the positive rows. The positive and negative lists are disjoint after cancellation. Apply SP1 with weights \((1,1,-1)\), all three original input caps and the original numerical shell to the remaining exact outputs. It counts at most \(O(n^{r/2+\eta}X^\varepsilon)\) completions. This yields FC2, summing over tail records, the one divisor \(F\), and the \(q-1\) anchored negative rows.

FP1 transfers FC2 to actual equal modular products of any admissible fixed word. With \(\mu\) small, \(R\) large and allowances small, the bound is \(n^{m/3-1/6+\delta}X^\varepsilon\) for every fixed \(\delta>0\) on eligible positive tuples, with choices depending on \(m,\delta\) and profile. Thus outside a power-small set of positive tuples, even competitors with different complete targets cannot give the same product, except permutations. The deletion loss can depend on \(m\) and profile. SWP1 already gives a different sparse-collision theorem for globally moving targets and arbitrary nonzero fixed signed weights. FC2 is confined to the equal unweighted products but gives a fixed limiting \(1/6\) gap for every fixed \(m\) on its generic-positive sector. FC1 is an unconditional cap-free bound on exact common-target energies. A separate upper bound on the occupied product image would be needed for a counting saving from these moments.

[Probe records and exact word lifting FP1](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#finite-probes-by-möbius-inversion)
· [Original mixed three-slot bound SP1](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#mixed-sign-completion-using-the-three-radical-caps)
· [Threefold coloured energy CE2](ABC_SAME_TARGET_COLOURED_PRODUCT_ENERGY.md#a-free-support-discount-for-common-target-threefold-energy)
· [Moving-target complementary energy SWP1](ABC_CROSS_TARGET_COMPLEMENTARY_WORDS.md)

## FC3. A nearly complete common numerical target lock

The common-target argument has a moving-target sector in which no third-slot completion theorem is needed. Assume that every positive target is divisible by one prescribed integer \(D\ge X^{1-\theta}\). Fix exact \(U,V\) on the signed lists and retain all original caps. A positive tuple is assumed to have at most \(\xi\) normalized free-support overlap per slot with the union of the other positive free-role supports. The condition does not count target overlaps. Suppose every candidate unanchored negative row supplies a cheap full divisor \(G\) of one of its original free roles, with \(\operatorname{rad}G\le n^\mu\) and \(G\ge X^\nu\). Take fixed choices satisfying
\[
 (m+1)\theta+R\xi<\omega<\nu/2.
\]
Then the nontrivial exact list-pair count is
\[
 O(T^{m-1}n^{\mu+2m/R}X^\varepsilon)
\]
where \(T\) bounds the size of the dictionary for freely fixed negative or canceled rows. The assumptions can also be imposed by prescribing the actual divisor catalogue, replacing \(\mu\) with its conditional price.

Cancel shared rows as in FC1, and fix \(q-1\) negatives if the remaining length is \(q\ge2\). Buy \(G\) on the one unanchored negative and fix every positive summand tail at raw cost \(2q/R\). For the \(A\) and \(B\) divisibility, replace \(\prod c_i\) in the exact clearing equations by \(\prod(c_i/D)\), since the marked positive summands are all coprime to \(D\). For each candidate positive row, the amount of both summands accounted for by this unknown target-cofactor product is at most \(q\theta\) in total numerical weight, as the candidate summands have disjoint supports. On private positive non-tail primes allocate available same-role negative factors as in FC2 after accounting for the known numerator and this positive-cofactor contribution. The logarithmic losses \(L_i\) charged to the unknown parts of the final negative summands satisfy \(\sum_iL_i\le\log_X(A'B'/G)\le2-\nu\), because allocation uses disjoint prime sets between distinct positive slots. At most one \(L_i\) exceeds \(1-\omega\). Adding the known \(D\) to the two known positive summand divisors gives a coprime triple of divisors whose product is at least
\[
 a_i b_i X^{\,\omega-(q+1)\theta-R\xi}
\]
on every good positive slot. This exceeds the two-dimensional height determinant on the dyadic positive summand boxes. Recover \(q-1\) positives by that determinant and divisor enumeration. The final positive and negative are completed by two affine ratio equations. This proves the claimed bound; \(q=1\) is immediate. A family-wide positive-support privacy assertion still has to be proved in any application conditioned on a lock \(D\), and the strength of the numerical lock is essential in the stated comparison.
