# Averaged full-modulus records and affine residue graphs

A nonzero cross determinant can pass only subpower many distinct full-modulus records, even when those records overlap. This gives a uniform averaged moment bound for actual signed root ratios and scalar ray tests. The gain from this averaging has an explicit limit once the record cost is included. A separate graph bound uses the three original radical caps to handle actual moving affine targets.

## A divisor bound for the moduli in a collision

Work at height \(X\), with \(n=X^L\) and the budget \(L>0\) fixed before the height limit. For two different positive reduced rational products \(P/Q,P'/Q'\) of polynomial height, agreement modulo a polynomial-height integer \(M\) at which both denominators are units implies
\[
                              M\mid PQ'-P'Q .
\]
Thus only \(X^{o_X(1)}\) distinct \(M\) can pass the test. The same holds for a vector of boundedly many prescribed moduli, all of which must pass. Repetitions or overlapping prime supports between the moduli do not affect this bound. At fixed anatomic depth an indexed full-divisor record has subpower multiplicity over its modulus, since each indexed factor divides that modulus. Arbitrary additional labels are separate records.

The corresponding assertion holds for the scalar ray tests in (PSR3). Fix the field \(K\) of bounded degree \(d\), the selection of degree-one prime ideals, and the power \(e\). For two different products \(\mathfrak A,\mathfrak B\) of polynomial ideal norm, a passing test gives a generator \(\gamma\) of \((\mathfrak A/\mathfrak B)^e\) whose centered logarithms are bounded by the width \(W\) of the test. There are only \(O_{d,W}(1)\) candidates for \(\gamma\) as the modulus changes, by the unit argument in (PSR3). Even \(W=O(\log X)\) gives only \((1+\log X)^{O_d(1)}\) candidates by uniform separation of bounded-degree unit log vectors. This follows from the bounded minimal polynomial coefficients of any ratio unit whose conjugates have bounded magnitude.

For one such \(\gamma\), clear its denominator as in (PSR3). The difference \(\Delta\) of two unequal conjugates of the resulting integer is nonzero because a rational generator would give \(\mathfrak A=\mathfrak B\) under the common orientation. Put \(D=[E:\mathbb Q]\), where \(E\) is a normal closure. Every full scalar modulus \(M\) passing with \(\gamma\) divides \(\lvert N_{E/\mathbb Q}\Delta\rvert\). This also holds for the omitted-degree-one-arm test in degree \(d\ge3\). For each \(p^a\parallel M\), the norm exponent is at least \(Da(1-2/d)\ge a\) by the fraction of places retained for both conjugates. The norm has polynomial height. Its divisor count and the subpower list of possible generators prove the scalar ray assertion. This uses no bound on the least common multiple of all tested moduli. [Projective and scalar ray collision tests (PSR1–PSR4)](ABC_PRIVATE_PACKET_RAY_SIEVE.md)

## Weighted moments over a modulus dictionary

Let \(R\) be a record with finite distribution \(\pi\). In applications this is the empirical distribution of a deterministic record on an actual family of primitive rows. For a group \(G\), let the tested projection of \(R\) be a bounded vector of actual moduli denoted by \(M_G(R)\), with
\[
                  \sup_s\Pr_\pi(M_G(R)=s)
                        \le n^{-\beta_G+o_X(1)}.               \tag{AM1}
\]
Extra residue targets need not be included in this projection for the moment bound, since the difference of two tests at the same record cancels the targets.

Suppose the entries in group \(G\) range over an ambient dictionary of size at most \(n^{a_G+o_X(1)}\), and over an arbitrary subdictionary at record \(R\) of size at most \(n^{c_G+o_X(1)}\). Take fixed positive integer moment order \(h_G\). For signed root ratios assume, as in (JD1), that identical collapsed products have subpower many ordered preimages once one product is fixed. Global coloring by exact root position ensures this, and the number of positions is fixed at the height limit. For ray ideal products the subpower bound follows from the ideal factorization count in fixed degree. For an off-diagonal pair assume the signed ratios or oriented ideals are formed from the same ambient entries before specifying the moduli. The preceding divisor argument then applies to a fixed off-diagonal pair. If several different tests are used in the projection, each counted modulus must divide a polynomial-height witness for the pair, as happens when the whole projection tests the same signed ratio or oriented ideal product. The field in a ray test is fixed across the records being averaged.

Let \(\mathcal E_G(R)\) be the \(2h_G\)-th character moment at this record. For a signed monomial take uniform character measure on the product of the unit groups; for a scalar ray test use the probability Fourier weights of the supported nonnegative kernels in (PSR4). The dictionary can have arbitrary record-dependent amplitudes of absolute value at most one. Expansion and absolute values give a nonnegative majorant whose diagonal part at every record is at most \(n^{h_Gc_G+o_X(1)}\). Each off-diagonal ambient pair can contribute only on subpower many distinct projections. Therefore
\[
       \mathbb E_\pi \mathcal E_G^{\rm off}(R)
                       \le n^{2h_Ga_G-\beta_G+o_X(1)},          \tag{AM2}
\]
where \(\mathcal E_G^{\rm off}\) denotes the off-diagonal part of that majorant. The same proof works if the whole norm packet dictionary changes arbitrarily inside its fixed ambient dictionary. In particular, for \(m\) distinct full modulus records given equal weight, an ambient dictionary of \(B\) entries has mean \(2h\)-moment
\[
                         X^{o_X(1)}(B^h+B^{2h}/m).
\]
Unlike (PSR4), this bound uses the number of distinct moduli and does not require their coprimality.

Now suppose the row count at a fixed \(R\) is bounded by the product of the corresponding \(2h_G\)-moment norms, with \(\sum_G(2h_G)^{-1}\le1\), up to an exponent allowance \(\epsilon_0 L\). This is the usual joint orthogonality and Hölder estimate for simultaneous unit congruences in (JD1). It applies to actual rows even with further coupled restrictions across groups, since one can enlarge to all congruent choices in the product of the local dictionaries. In a count with amplitudes, this hypothesis is imposed on the actual row majorant; the arbitrary amplitudes in (AM2) alone do not give such a majorant. For available ordinary diagonal estimates use \(c_G/2\) as the corresponding cost.

If the support of \(R\) has size at most \(n^{b+o_X(1)}\), then (AM2) and Markov yield the count
\[
                   N \ \le\
               n^{\,b+\sum_G
                   \max\{c_G/2,\ a_G-\beta_G/(2h_G)\}
                        +\epsilon_0+\eta+o_X(1)}               \tag{AM3}
\]
for every fixed \(\eta>0\), with the stated substitution on ordinary diagonal groups. In fact all but a power-smaller fraction of the record law satisfy the pointwise off-diagonal majorants with \(n^\eta\) allowance, adjusted before the inner limit for the fixed number of groups. Apply the per-record Hölder bound there and sum over the support. At fixed truncation the exceptional fraction tends to zero, so this bounds the size of the actual family when \(\pi\) is its empirical record law. The local trivial cap \(c_G\) can also be used if it is smaller than the displayed maximum. [Signed joint dictionary moments (JD1)](ABC_JOINT_FREE_LETTER_DICTIONARIES.md#a-signed-joint-dictionary-estimate)

One can bin the realized record surprises \(-\log_n\pi(R)\) to use an empirical support cap. A bin with surprise \(b+O(\eta)\) has at most \(n^{b+O(\eta)}\) possible records. If it retains a subpower fraction of the rows, passage to its conditional law preserves the corresponding min-entropy estimates up to the same allowances. The ambient dictionaries in (AM2) must be chosen uniformly within the bin.

## The cost of averaging and the odd-kernel moment

Since \(M_G(R)\) is a projection of a record with support at most \(n^{b+o_X(1)}\), we have \(\beta_G\le b\). On the off-diagonal branches the total exponent reduction from (AM1) is therefore at most
\[
                      \sum_G\beta_G/(2h_G)\le b .
\]
Consequently that reduction by itself cannot exceed the recording cost. In (AM3) a strict improvement requires an additional restriction, such as a genuine diagonal bound or a compressed ambient dictionary.

For example, at balanced caps suppose all free groups other than one have diagonal moments, and the remaining group of raw cost \(f\) has ambient and local cost \(f-s\). Use \(h_G=1\) on this group and suppose the tested modulus record has both support exponent and min-entropy exponent \(b\), up to the inner allowances. Equation (AM3) gives the normalized count coefficient
\[
               \tfrac13+\max\{b-s/2,\ (b+f-2s)/2\}.            \tag{AM4}
\]
This is strict below \(1/3\) when \(b<s/2\) and \(b+f<2s\), with fixed margins. For the moving square kernels compressed by (JD5), \(s=\theta/2\), so the second condition would require \(b+f<\theta\). A nonempty high-torsion sector at this cap has \(\theta\le f+o(1)\) by the elementary class-number bound \(\lvert\mathrm{Cl}(K)\rvert\ll_\epsilon |\operatorname{Disc}(K)|^{1/2+\epsilon}\) on the two quadratic factors. Both discriminants are \(O(n^f)\). Thus first-moment compression and the arithmetic mean (AM2) alone cannot eliminate the residual high-torsion sector of (JD6). [Moving square kernels and the unresolved moment (JD5–JD6)](ABC_JOINT_FREE_LETTER_DICTIONARIES.md#the-unrestricted-kernel-group-in-this-mixed-signature)

## An actual affine residue graph

Fix a positive integer \(k\ge2\), constants \(B,H,w,S>0\), and a radical cost cap \(\mu\ge0\) with
\[
                          kw>H,\qquad k\mu<1/6 .
\]
Take \(m\) pairwise nonproportional forms \(T_i(x)=a_i x+b_i\) with positive integer coefficients at most \(X^B\), and \(N\) distinct positive integer queries in an interval of width at most \(X^H\). A graph edge \((i,x)\) requires \(a_ix+b_i=T_i(x)\le X\) to be primitive and each of its three complete values to have radical at most \(n^{1/3+\delta_0}\). Also require a full divisor of \(T_i(x)\) with height at least \(X^w\) and radical at most \(n^\mu\), whose prime exponents satisfy \(L v_p(T_i(x))\le S\) on its support. One may choose this witness separately for each edge.

There are \(\delta_0,\delta,L_0>0\) such that for \(0<L<L_0\) fixed the number of edges satisfies
\[
                |\mathcal E|
        \ll_k N+
          m N^{1-1/k}\,
           n^{\max(k\mu,\,1/6-\delta)/k+o_X(1)} .             \tag{AM5}
\]
The inner allowance can be chosen arbitrarily small. Uniformity is in \(m\), the forms, their libraries of witnesses, and the query set, at fixed parameters.

Choose \(0<\xi<1/12\) small enough that \(k(w-S\xi)>H\). For each fixed \(k\)-tuple \(I\) of distinct anchors put all primes of its slopes and its nonzero pair determinants into a pool \(\mathcal P_I\) of polynomial-height prime product. For any \(i\in I\), consider its neighbors \(x\) for which \(T_i(x)\) has radical mass at least \(\xi\) in \(\mathcal P_I\). Their target radical outside this pool is at most \(n^{1/3+\delta_0-\xi}\). The witness forces \(T_i(x)\ge X^w\), so each dyadic shell \(Y/2<T_i(x)\le Y\) has \(\log_X Y\ge w\). Cover the shell exponents by finitely many fixed intervals \([h_-,h_+]\subseteq[w/2,2]\) with \(h_+/h_-\) sufficiently close to one in terms of \(\xi\). Use the fixed budget \(L'=L/h_-\) in each interval. The raw caps and outside cap in base \(X^L\) also hold in base \(Y^{L'}\), since \(Y^{L'}\ge X^L\). Apply (FL1) on that shell with \(b_i\) fixed, the target pool \(\mathcal P_I\), and the other two raw caps. The pool product is still polynomial in \(Y\), uniformly in \(I\). Take \(L_0\) small enough for all the \(L'\) and \(\delta_0\) and allowances sufficiently small in terms of \(\xi\). Replacing \(Y^{L'}\) by its upper bound \(n^{h_+/h_-}\) and summing \(O(\log X)\) shells gives at most \(n^{1/6-\delta+o_X(1)}\) such neighbors, for one \(\delta>0\) independent of \(I,i\). [Fixed-letter count with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

For every other common query, strip all primes of \(\mathcal P_I\) from the chosen witnesses. This loses at most \(S\xi\) in numerical height on each. The remaining moduli are pairwise coprime, since a prime dividing two corresponding form values divides the determinant of those forms. Each slope is a unit at the remaining modulus, so the congruence \(T_i(x)\equiv0\bmod M_i\) fixes \(x\) to a single residue class. For fixed stripped moduli their product exceeds \(X^H\), giving at most one query by the Chinese remainder theorem. Each stripped modulus has at most \(n^{\mu+o_X(1)}\) possibilities by its radical cap and its polynomial height. Thus every \(k\)-tuple has at most \(n^{k\mu+o_X(1)}+k n^{1/6-\delta+o_X(1)}\) common queries. Sum over anchor tuples, then use the power mean for the query degrees at least \(2k\) to obtain (AM5). The other query degrees contribute at most \(2kN\). Each pool is used only for its own fixed-size tuple. The aggregate product of all anchor determinants is unrestricted.

For example, with \(N=n^{1/6+o_X(1)}\) and power-many anchors, (AM5) gives a fixed power saving in the graph density. With polylogarithmically many anchors it bounds the average query degree by a constant. A bounded compact valuation band strictly above scaled degree \(6H\) supplies the witness slope for this application if it contains enough radical mass in sufficiently small indivisible packets. Arbitrary balanced profiles need not supply that band. The inequality describes the additional height and cost criterion in an actual graph with all three radical caps.


### Projective query pairs

There is a projective version with \(1/3\) in place of
\(1/6\) in the hypothesis \(k\mu<1/6\) and the second term
of (AM5). Use \(N\) distinct positive primitive pairs \((u,v)\)
whose pairwise nonzero cross determinants have absolute value
at most \(X^H\). A graph edge now gives the actual primitive
triple
\[
                  a_i u+b_i v=T_i(u,v)\le X
\]
with the same three radical caps and full-divisor witness
conditions. The coefficients and the pair determinants of the
forms have the same polynomial height bounds.

For a fixed \(k\)-tuple of anchors, put the primes of both
coefficients of each form and the nonzero pair determinants in
the pool. Apply the effective radical count (EB4) to the
original triples for one fixed anchor. Its total effective
cap drops from \(1+3\delta_0\) to at most \(1+3\delta_0-\xi\)
on the bad branch. The same shell mesh as for (AM5) gives
at most \(n^{1/3-\delta+o_X(1)}\) bad queries. The map
from \((u,v)\) to the two complete summands is injective
for a fixed anchor. For a good query, strip the witnesses
as before. If a prime divides two form values on a
primitive pair, it divides the corresponding determinant.
Thus the retained moduli are pairwise coprime. At each
modulus the two coefficients are units, and the form fixes
the primitive projective class. Distinct pairs with the
same classes have a nonzero cross determinant divisible
by the product of these moduli, which is greater than
\(X^H\). Therefore each set of modulus records fixes
at most one query. The same codegree argument proves
\[
               |\mathcal E|
         \ll_k N+mN^{1-1/k}
                     n^{\max(k\mu,\,1/3-\delta)/k+o_X(1)}
\]
and gives a strict density saving at \(N=n^{1/3+o_X(1)}\) for
power-many anchors. The projective height condition may
come from an angular interval. No estimate on all projections
is required.
[Effective radical count (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
· [Projective angular comparison (PSR1)](ABC_PRIVATE_PACKET_RAY_SIEVE.md#projective-tests-in-an-angular-interval)
