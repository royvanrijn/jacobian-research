# Prime membership and projective incidence after least repair

Repairing all three letters separately sends a primitive row to a projective point. Within a row the three repaired letters are still pairwise coprime. In a polynomial relation across several such rows, a prime can enter at most one coordinate of each row. Restricting the permissible prime membership patterns before computing the invisible directions gives stronger incidence bounds. The same argument works for rational monomials of the repaired letters, including coefficients attached to other rows in a fixed motif.

## Counting with admissible prime membership

Fix the number of rows \(h\), a bound \(B>0\), and fixed lists of monomial supports. Let \(n=X^L\), with the height limit taken at fixed \(L>0\). Take \(K>0\) and \(M=\lceil K/L\rceil\). The input rows are positive primitive \(a_i+b_i=c_i\le X\) with individual radicals at most \(n^B\). Write
\[
              H_{i,j}=H_M(A_{i,j}),\qquad
              A_{i,1}=a_i,\quad A_{i,2}=b_i,\quad A_{i,3}=c_i .
\]
Only the coordinates used below need their radical caps. For fixed vectors \(\mathbf w_{i,j}\in\mathbb Z^q\), use \(q\) rational observables
\[
                  t_\ell=\prod_{i,j} H_{i,j}^{w_{i,j,\ell}}
                       \qquad(1\le \ell\le q).
\]
Call a base coordinate \((i,j)\) active if \(\mathbf w_{i,j}\ne0\), and let \(m_0\) count the active coordinates. The possible nonzero membership vectors are
\[
 {\cal A}
 =\left\{\sum_i\mathbf w_{i,j_i}: j_i\in\{0,1,2,3\},
                               \ \mathbf w_{i,0}=0\right\}
                        \setminus\{0\}.
\]
The choice of one index per row follows from primitivity. The entries of a vector in \({\cal A}\) can be signed or can have magnitude greater than one.

Prescribe any fixed number of nonzero Laurent polynomials \(G_b(\boldsymbol t)\), each a sum of at least two distinct monomials with nonzero rational coefficients. Require each equality \(G_b(\boldsymbol t)=0\) to be nondegenerate as an identity of its displayed terms. A membership vector \(\boldsymbol v\) is invisible to one \(G_b\) if \(\boldsymbol e\cdot\boldsymbol v\) is constant over all exponent vectors \(\boldsymbol e\) of that polynomial. Put
\[
 V=\operatorname{span}_{\mathbb Q}
       \{\boldsymbol v\in{\cal A}:
                        \boldsymbol v\text{ is invisible to every }G_b\}.
\]
For every \(\eta>0\), at all sufficiently small fixed \(L>0\), the observable tuples lie in at most
\[
                         O(n^{m_0/K+\eta}X^\varepsilon)      \tag{HM1}
\]
classes of the relation \(v_p(\boldsymbol t/\boldsymbol t')\in V\) for every prime \(p\). Here \(v_p\) is applied coordinatewise. The exponents and vector lists are fixed before the limits, and the constant is uniform in the coefficients. The same proof allows any smaller exponent for the joint tail dictionary when it is available.

To prove this, let \(W_{i,j}=\prod_{p\mid A_{i,j}}p^{(v_p(A_{i,j})-M)_+}\) for each active coordinate. Its radical is at most \(X^{1/(M+1)}\). Fix all these tails at cost at most \(X^{m_0/(M+1)+\varepsilon}\) by radical enumeration and the uniform fixed-support count. For each \(\boldsymbol v\in{\cal A}\) form \(x_{\boldsymbol v}\) from the primes having this membership vector. These are pairwise coprime squarefree integers of size at most \(n^{Bm_0}\). The primes with vector zero contribute only through the already fixed tails. The observables have the form
\[
                   t_\ell=w_\ell\prod_{\boldsymbol v\in{\cal A}}
                                      x_{\boldsymbol v}^{M v_\ell}
\]
with fixed rational \(w_\ell\). Each equation factors off the variables invisible to it. Multiplying by a further monomial clears negative exponents. Each remaining variable has spread at least \(M\), and the total number of variables is bounded in the fixed data.

For at least three displayed terms with visible variables, use the coefficient-uniform nondegenerate power bound (EB1). It permits terms with the same resulting exponent vector and proves the requisite saving for the visible variables. Equivalently, group terms with the same exponent vector. With a visible variable there are at least two groups, and nondegeneracy forbids a group coefficient from vanishing. If only two monomials are left, take their ratio. The resulting rational constant has polynomial height whenever a solution occurs, because it is a product of the actual visible squarefree integers raised to fixed powers at the fixed budget. Each prime of that constant can be assigned to at most \(|{\cal A}|\) slots, and every prime of a visible slot must occur in the constant. This gives \(X^\varepsilon\) choices uniformly even when the original coefficients are unbounded. For a single visible variable the univariate bound also applies. An equation with no visible variables costs nothing to record. The resulting exponent for the visible variables in each equation is at most \(\eta\) divided by the number of equations, in base \(n\), at sufficiently small fixed \(L\). Fix the union of these visible records. A variable still free is invisible to every equation, so the remaining differences of prime valuation vectors lie in \(V\). This proves (HM1). [Power release (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md) · [Unrestricted membership (HP1)](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md#prime-packets-invisible-to-a-polynomial-relation)

A set of coordinate projections separating \(V\) may be used to count whole observable vectors inside a class. If fixed observed coordinates give prescribed repairs of two different letters in one primitive row, their supports determine those two original integers with \(X^\varepsilon\) choices, and thus the row with \(X^\varepsilon\) choices. The same count applies if one known observable is a monomial of at least two repairs from a single row, all its exponents on those letters nonzero. The numerator and denominator support determines the union of the supported primes because the original letters have disjoint supports. One can enumerate the original letters on that union directly by the fixed-support bound. No lower numerical height is required.

## Rational observables indexing entire primitive rows

Choose for each slot \(i\) a fixed integer triple \(\boldsymbol d_i=(d_{i,1},d_{i,2},d_{i,3})\) with at least two nonzero entries, and put
\[
                         t_i=\prod_j H_{i,j}^{d_{i,j}} .
\]
For example \(H_M(a_i b_i)=H_{i,1}H_{i,2}\) and \(H_M(a_i b_i)/H_M(c_i)^2\) both qualify. If each input row comes from a family of \(T\) primitive rows, then the number of row tuples satisfying the prescribed nondegenerate relations is
\[
                            O(T^r n^{m_0/K+\eta}X^\varepsilon),
                           \qquad r=\dim V .                \tag{HM2}
\]
Indeed select \(r\) scalar slots whose projection separates \(V\), and prescribe original rows in them. The class and these observables determine the others. Each remaining scalar indexes its entire row at subpower cost. When each observable is the repair of a pair product, \(m_0=2h\). In this case the gain from a dimension drop is a full row exponent, without a complete-target fibre factor.

In a homogeneous linear relation with these \(h\) scalar slots and all coefficients nonzero, invisibility means all coordinates of \(\boldsymbol v\) agree. Thus \(r=1\) if the nonzero weight sets \(\{d_{i,j}:d_{i,j}\ne0\}\) have a common element and \(r=0\) otherwise. A specified nonzero shift makes \(r=0\) in either case. If the rational weights have polynomial height and the denominator primes lie in a fixed polynomial-height pool, a varying nonzero shift with outside numerator radical at most \(n^b\) costs at most \(n^{b}X^\varepsilon\). In fact each observable has polynomial numerator and denominator, and the shift denominator is supported on the primes of the prescribed weights and the negative-weight base coordinates. To make this denominator support fixed, either require the observables to be integral, or prescribe a polynomial-height pool containing all their denominator primes. The shifted count with that requirement applies to a nondegenerate sum. If the weights and observables are positive, nondegeneracy is automatic.

For a single positive observable type, the dictionary has size \(T X^{-\varepsilon}\) or larger, up to the implied constant. With a pair-product type it has the same ordered equal-sum estimate as (HP4) with tail charge \(4h/K\) for two \(h\)-term sums. For \(T\ge n^d\) with fixed \(d>0\) and \(K\) sufficiently large, the \(h\)-fold sumset of that dictionary has cardinality asymptotic to the dictionary size to power \(h\) divided by \(h!\). This counts values indexing two original letters with no need to estimate the fibre of one target. [Energy argument (HP4)](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md#shifted-sums-and-true-additive-diagonals)

## Bounded-degree curves in the projective repair image

For one primitive row all permissible membership vectors of the repaired coordinate triple are singletons. Let \(G(U,V,W)\) be any nonzero homogeneous polynomial of bounded degree with rational coefficients. Then
\[
 \#\{(a,b,c):G(H_M(a),H_M(b),H_M(c))=0\}
                         \ll n^{3/K+\eta}X^\varepsilon.     \tag{HM3}
\]
Partition the polynomial terms of a solution into minimal vanishing blocks. They all have at least two terms, since the point has positive coordinates. Each block varies in at least two coordinates because it is homogeneous and its terms have distinct exponent vectors. For any one such block, (HM1) fixes the repairs of the varying coordinates once the visible records are fixed. This uses the singleton membership list, and its cost is bounded by the right side of (HM3). Two known repairs recover the row at subpower cost. Only finitely many blocks are possible. The bound is uniform over all polynomials of the given degree bound by taking a finite union of allowable monomial supports. The same affine statement holds on the complement of the positive rational coordinate slices that are factors of \(G\). In fact if the visible-coordinate union of a vanishing partition has size at most one, every block, after removing a monomial, is a polynomial in that same coordinate. It vanishes at the slice value, so \(G\) vanishes identically on the corresponding slice.

Consequently any family of bounded-degree curves with a dictionary of at most \(n^\gamma\) coefficient records covers at most \(n^{\gamma+3/K+\eta}X^\varepsilon\) rows in the projective repair image. This applies to a moving projective change of coordinates chosen from that many maps, when the images are constrained to a fixed curve. The projective points from the original integer triples all lie on the additive line. The actual separate repairs of a near-saturating family cannot all be carried from that line by a single projective linear map.

## Collinear triples of distinct repair vectors

Let \(P_i=(H_{i,1},H_{i,2},H_{i,3})\) for three original rows. Among ordered triples of input rows with pairwise distinct repair vectors, one has
\[
               \#\{\det(P_1,P_2,P_3)=0\}
                          \ll T n^{9/K+\eta}X^\varepsilon . \tag{HM4}
\]
Each repair vector is primitive since its coordinates are pairwise coprime, so distinct vectors are distinct projective points. In particular, for a family of exponent bounded below by a fixed positive number, there is a large subset of the projective repair dictionary with no collinear triple. If that dictionary has size \(D\), then \(T\ll D X^\varepsilon\) by the uniform support count. Independent sampling and deleting one point for each surviving triple gives a subset of size \(\gg D n^{-9/(2K)-\eta}X^{-\varepsilon}\), after resizing the allowances. Choose sampling probability a sufficiently small constant times the inverse square root of the per-vertex bound \(n^{9/K+\eta}X^\varepsilon\).

Here is the determinant count, including its cancellations. A membership pattern selects at most one of the nine cells in each row. The six permutation terms of the determinant hit a total of twice the number of selected cells. If all six terms have the same hit count, the pattern is empty or selects three cells. In the latter case they must all be in the same column, since any pattern with two cells in different columns admits a permutation hitting both. Conversely a full column hits each term once. The invisible space is therefore spanned by the three column indicators and is separated by any one complete row. For a nondegenerate determinant relation, fixing one original row and applying (HM1) gives the right side of (HM4). The remaining repaired vectors fix the remaining primitive rows at subpower cost. The restriction on patterns matters here. In the unrestricted pattern space the invisible row and column indicators already span a space of dimension five.

For a degenerate determinant, first suppose some two terms cancel. Their permutations have opposite signs, so one is obtained from the other by swapping two columns on two rows. Dividing the common third cell gives equality of the cross products of the selected two entries on these rows. The entries of each pair on a row are coprime, hence the same two repaired entries occur on both rows. If their third entries are also equal, the full repair vectors coincide, excluded in (HM4). Otherwise subtract those two rows in the determinant. Its zero value forces the third row to have the same ratio of the selected two entries, hence the same two entries by coprimality. One original row recovers both others at subpower cost. There are only finitely many choices of the cancelled terms.

If no two terms cancel and the determinant is degenerate, a minimal vanishing partition has two blocks of size three. A pattern invisible to both blocks has twice its size divisible by three, so it is empty or selects exactly one cell per row. For three selected cells in distinct columns the six hit counts are \(3,1,1,1,0,0\), by counting fixed points of the relative permutation. For two selected cells in one column and the third in a different column, they are \(2,2,1,1,0,0\), as seen by moving the selected columns to \((1,1,2)\). Neither multiset can be constant on two blocks of size three. The simultaneously invisible patterns are again the full columns and the empty pattern. Apply (HM1) to the two nondegenerate block relations together. This proves (HM4). Allowing identical repair vectors adds \(\ll T^2X^\varepsilon\) tuples.

The same membership argument can eliminate all invisible patterns when the monomials of an equation have disjoint sets of individual repaired coordinates, all observed coordinates occur, and the number of terms exceeds the number of input rows. Every nonempty admissible pattern then meets at least one but misses another term. For instance a relation
\[
                    H_{1,1}H_{2,2}+H_{1,2}H_{2,3}
                                      =H_{1,3}H_{2,1}
\]
has at most \(O(n^{6/K+\eta}X^\varepsilon)\) ordered row pairs. The three-term relation is automatically nondegenerate.

## Numerator counts without a fixed denominator

Let \(P(\boldsymbol t)=\sum_\nu A_\nu\boldsymbol t^{\boldsymbol e_\nu}\) have fixed Laurent support with at least two terms. Assume the reduced rational coefficients have polynomial height in \(X\), uniformly. For an admissible direction \(\boldsymbol v\) define
\[
             \mu(\boldsymbol v)
                   =\min\bigl(0,\min_\nu \boldsymbol e_\nu\cdot
                                                       \boldsymbol v\bigr).
\]
Suppose whenever \(\mu(\boldsymbol v)<0\), exactly one term attains that minimum. Write \(P(\boldsymbol t)=u/z\ne0\) in lowest terms with \(z>0\). Assume the augmented relation between the terms of \(P\) and \(-u/z\) is nondegenerate. For an outside radical cap
\[
                    \operatorname{rad}_{\overline{\mathcal P}}(u)
                                      \le n^b
\]
with prescribed polynomial-height prime pool \(\mathcal P\), the number of observable valuation classes is
\[
                            O(n^{b+m_0/K+\eta}X^\varepsilon). \tag{HM5}
\]
Here the equivalence space is spanned by those admissible \(\boldsymbol v\) satisfying \(\boldsymbol e_\nu\cdot\boldsymbol v=0\) for every term. Consequently when this space is zero, the bound applies to whole observable tuples. For the row-indexing scalar observables of (HM2), a space of dimension \(r\) costs at most \(T^r\).

Fix the tails and enlarge the prime pool by their supports and by the numerator and denominator primes of the coefficients. The enlarged pool has polynomial-height product. For each actual observable and the denominator \(z\), fix its part supported on this pool. These are signed-valuation records of polynomial height and cost \(X^\varepsilon\) uniformly by fixed-support enumeration. Fix also \(u\), at cost \(n^bX^\varepsilon\). Form the packet roots \(x_{\boldsymbol v}\) using the primes outside the enlarged pool. At those primes term valuations are exactly \(M\boldsymbol e_\nu\cdot\boldsymbol v\). If their minimum is negative, the unique least valuation gives precisely \(M\mu(\boldsymbol v)\) for the valuation of \(P(\boldsymbol t)\); if it is nonnegative these primes do not enter its denominator. Therefore
\[
 \begin{split}
 t_\ell&=w_\ell\prod_{\boldsymbol v}x_{\boldsymbol v}^{Mv_\ell},
 &z&=z_0\prod_{\boldsymbol v}x_{\boldsymbol v}^{-M\mu(\boldsymbol v)},\\
 u&=\sum_\nu z_0 A_\nu\boldsymbol w^{\boldsymbol e_\nu}
           \prod_{\boldsymbol v}
              x_{\boldsymbol v}^{M(\boldsymbol e_\nu\cdot
                                            \boldsymbol v-\mu(\boldsymbol v))}.
 \end{split}
\]
Here \(w_\ell,z_0\) are fixed by the supported records. Multiplication by \(z\) preserves nondegeneracy. The nonzero constant term makes a packet invisible only if all its exponents on the right vanish. The unique-minimum hypothesis then forces all original term exponents to be zero. Apply (EB1) to the remaining visible packets, with the two-term and other terminal cases handled as in (HM1). This proves (HM5). The pool-supported parts can be fixed at subpower cost despite the moving denominator. The outside part of the denominator has been recovered as a monomial of the same variables.

In particular, take a linear sum of the row-indexing observables in (HM2) with nonzero weights. If their negative weights \(d_{i,j}\) have no common value across two different rows, the unique-minimum hypothesis holds. The invisible space in (HM5) is zero. For example one may use \(t_i=H_M(a_i)H_M(b_i)/H_M(c_i)^{j_i}\) with distinct positive integers \(j_i\). With positive sum weights the relation is automatically nondegenerate. If a family has \(T\ge n^d\) rows for fixed \(d>0\), all but a power-smaller fraction of its \(h\)-tuples have reduced numerator of outside radical greater than \(n^{hd-3h/K-\rho}\), for any fixed \(\rho>0\), after choosing the inner allowances below \(\rho\). When all tails vanish one can omit the \(3h/K\) in this threshold. This conclusion needs no restriction on the denominator primes of the sum.

Allowing repeated negative weights can invalidate such a count if the denominator cancellation is left uncontrolled. Use the same weight type \(H_M(a)/H_M(c)\) in two slots, take the same primitive row \(1+(R^M-1)=R^M\) in both, and let \(R\) range over odd squarefree integers up to \(X^{1/M}\). The two active tails vanish and the reduced numerator of the sum is always \(2\). The active radical caps hold for \(B\ge1/K\), yet there are \(X^{1/M+o(1)}\) such pairs. This tests the count with only the active caps, without a separate cap on the middle letter.

### Repeated observables and fresh prime support

There is a version allowing repeated negative weights when the shared input radical is small. Let \(J\) be the product of primes on at least two input rows. On tuples with \(J\le n^\rho\), fix it at cost at most \(n^\rho\) and add it to the recorded pool. For a linear sum of row-indexing scalar observables, a remaining packet is private to one input row. Any negative term valuation is therefore the unique minimum. The proof of (HM5) applies with just the private directions, even if the same rational observable is used in every slot. It gives exponent \(b+m_0/K+\rho+\eta\) for the count with a numerator radical cap \(n^b\) outside the original pool and \(J\), for a nondegenerate sum or affine shift. A constant term in the output expression has exponent zero and does not affect uniqueness of a negative minimum. For each nonzero private direction the augmented expression still has a nonzero exponent and the observed equivalence space is zero.

For example, use \(h\ge2\) rows in a dyadic shell, with total radical at most \(n\) per row, and set
\[
                        t=\sum_{i=1}^{h}\frac1{H_M(a_i b_i c_i)}
                                      =\frac{N}{D}
\]
in lowest terms. This is a positive sum and \(t\le2h/X<1\) for large \(X\). Both augmented relations with outputs \(t\) and \(t-1\) are nondegenerate. For \(t-1\), the only negative displayed term in the zero relation is the shift \(-1\). Each private input prime has a unique negative valuation in the sum, so remains in \(D\) and divides neither \(N\) nor \(D-N\). Consequently the numerator cap outside \(J\) and a prescribed pool is exactly the fresh cap outside all input primes and that pool.

For each fixed \(\tau>0\), take \(K\) sufficiently large in terms of \(h,\tau\), choose a sufficiently small \(\delta>0\), and take the height limit for sufficiently small fixed \(L\). If the input family has at least \(n^{1/3-\delta}\) rows, then outside a power-smaller fraction of its ordered \(h\)-tuples,
\[
 \operatorname{rad}(D)\ge n^{h-\tau},\qquad
 \operatorname{rad}_{(\mathcal I\cup\mathcal P)^c}(N),
 \ \operatorname{rad}_{(\mathcal I\cup\mathcal P)^c}(D-N)
                                        \ge n^{h/3-\tau},   \tag{HM6}
\]
where \(\mathcal I\) is the set of primes on any input row. Indeed, for each sufficiently small fixed \(\rho>0\), the prescribed-pool discount excludes, at a power saving, tuples with \(J>n^\rho\). Fix one row of an offending pair and discount its support on the other row. The reduced-budget bound also excludes rows of total radical below \(n^{1-\rho}\). These are the shared-prime and small-budget exclusions used in (NR2). The retained private primes have total normalized mass at least \(h(1-2\rho)\). For either numerator, the preceding fixed-\(J\) argument with \(b=h/3-\tau\) bounds the bad tuples by \(O(n^{h/3-\tau+3h/K+\rho+\eta}X^\varepsilon)\). Choose \(\rho,\eta,3h/K,h\delta,\varepsilon/L\) with total below a fixed fraction of \(\tau\), and choose \(\delta\) also below the positive margins in the exclusions. Each \(t_i\) records the support of its three input letters, hence lifts at subpower cost. The output is the positive primitive additive triple \(N+(D-N)=D\), and its three radicals are pairwise coprime. [Private and shared input primes (NR2)](ABC_PROJECTIVE_SUM_NUMERATOR_RADICALS.md#private-denominator-primes-and-support-records) · [Prescribed-pool discount (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The incidence results count primitive quotients, with a fixed least repair attached to each chosen letter. Their arithmetic input is the coefficient-uniform power release, together with the restriction on where each prime can occur. The external generalized Fermat result works with an additive relation on the original full numbers. The bounds here apply to arbitrary bounded-degree relations on separately repaired numbers and allow nontrivial invisible directions across rows. They leave open whether an actual balanced family forces enough of the latter relations or cheaply supported shifted outputs. [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1) · [Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)