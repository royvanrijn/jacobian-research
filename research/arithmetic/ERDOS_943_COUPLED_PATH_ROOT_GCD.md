# Coupled path-root gcds in actual rectangular six-cycles

Let a simple bipartite graph have squarefree positive left kernels b<=K1, right kernels d<=K2, and an edge for a positive representation
\[
 \kappa_1 a^2b^3+\kappa_2 c^2d^3=n,\qquad
 \gcd(\kappa_1a^2b^3,\kappa_2c^2d^3)=1,\quad a\le R_1,\ c\le R_2. \tag{1}
\]
Select ONE actual representation per edge. All labels below refer to that selection. Coefficients are positive and fixed in a graph; constants will be uniform in them. Assume n>=2, all four caps are between 1 and n, and
\[
                         2(R_1R_2)^2<n.                       \tag{2}
\]
Only simple six-cycles are counted. Any usual labelled/unlabelled convention changes an absolute factor.

**Theorem 1 (Uniform coupled-path cycle tail).** Under (1)--(2), for every Y>=1 the number of actual selected simple six-cycles with alpha^4 gamma>=Y at some shared-right path is
\[
 O_\varepsilon(n^\varepsilon K_1R_1^2R_2^2/Y).
\]
Moreover the sum of alpha^4 gamma over all shared-right paths in all such cycles is O_epsilon(n^epsilon K1 R1^2 R2^2). Here alpha and gamma are the two gcds in (3) at the SAME path. The constants are uniform for each fixed positive coefficient pair. Swapping summands gives the shared-left statement. Sections 1--3 prove the theorem.

## 1. A lifted bilinear constraint

At a shared-right path with distinct b1,b2 and common d, set
\[
 a_i=\alpha r_i,\ c_i=\gamma s_i,\quad
 \alpha=\gcd(a_1,a_2),\quad\gamma=\gcd(c_1,c_2).                \tag{3}
\]
Every root, kernel and coefficient on an edge is a unit modulo n. Within-edge coprimality also gives gcd(a_i,c_i)=1, hence gcd(alpha,gamma)=1, gcd(r1,r2)=gcd(s1,s2)=1, and gcd(alpha,kappa2*d)=1. Actual subtraction gives a nonzero integer
\[
 \ell=(c_1^2-c_2^2)/\alpha^2.
\]
If c1=c2, squarefree prime-valuation parity in a1^2 b1^3=a2^2 b2^3 would force b1=b2, which is excluded. Since alpha and gamma are coprime,
\[
 \ell=\gamma^2\ell_0,\quad
 \alpha^2\ell_0=(s_1-s_2)(s_1+s_2),\quad \ell_0\in\mathbb Z\setminus\{0\}.
                                                                    \tag{4}
\]
Thus
\[
 \alpha\gamma<R_2,\qquad \gcd(s_1-s_2,s_1+s_2)\mid2.           \tag{5}
\]
Indeed alpha^2<=|s1^2-s2^2|<(R2/gamma)^2. This is an integer factorization, stronger than a congruence. In particular an odd prime power p^e|alpha obeys p^(2e)<2R2/gamma: its square divides one of the two factors in (4). If 2^e|alpha, the analogous bound is 2^(2e)<4R2/gamma, since s1,s2 are odd and one of their sum/difference has valuation one.

## 2. Actual codegrees and recovery

For fixed distinct b1,b2, elimination gives
\[
 b_1^3(a_1c_2)^2\equiv b_2^3(a_2c_1)^2\pmod n.
\]
The reduced positive fraction u/v=(a1*c2)/(a2*c1) has unit numerator/denominator <=R1*R2. It lies in at most 2^(omega(n)+1) unit square-root classes. A unit has at most two roots modulo an odd prime power and four at a power of two. Two distinct such fractions in one class have nonzero cross difference divisible by n, but of absolute value <=2(R1R2)^2<n. Thus at most one fraction occurs per class.

In (3), u=r1*s2 and v=r2*s1 are already coprime: the within-edge gcds gcd(r_i,s_i)=1 and pairwise gcds of r1,r2 and s1,s2 exclude every possible common prime. Divisor-many factorizations recover r1,r2,s1,s2 from u,v. The identity
\[
 \kappa_1\alpha^2(b_1^3u^2-b_2^3v^2)=n(s_2^2-s_1^2)
\]
determines alpha at most once; its bracket is nonzero by squarefree parity. The first edge then determines
\[
 \gamma^2d^3=
 (n-\kappa_1b_1^3\alpha^2r_1^2)/(\kappa_2s_1^2).
\]
A positive integer has at most one expression as a square times a squarefree cube. Invalid candidates are discarded. All integers factored are polynomially bounded in n. Therefore either same-side pair has O_epsilon(n^epsilon) actual common neighbors, uniformly in coefficients.

A specified selected shared-right path has O_epsilon(K1*n^epsilon) simple-cycle continuations: choose the third left kernel, then use two codegrees. A fixed kernel triple (b1,b2,d) specifies at most one selected path.

## 3. Coupled coefficient frequency

Fix real A,C>=1. Paths with alpha>=A and gamma>=C are empty if AC>=R2. Otherwise
\[
 r_i\le R_1/A,\quad 0<|\ell|\le L=R_2^2/A^2,
\]
and gamma^2|ell with gamma>=C. The number of possible nonzero ell is at most
\[
 2\sum_{C\le h\le\sqrt L}\lfloor L/h^2\rfloor\ll L/C.         \tag{6}
\]
For C>=2 use an integral comparison; for 1<=C<2 use convergence of sum h^-2. An empty range contributes zero. Overcounting ell with several square divisors is harmless.

The actual primitive kernel triple solves
\[
 \kappa_1r_1^2 b_1^3-\kappa_1r_2^2 b_2^3+\kappa_2\ell d^3=0. \tag{7}
\]
All coefficients are nonzero. Each b_i is coprime to d, so gcd(b1,b2,d)=1. Heath-Brown, arXiv:2601.07817v1, Theorem 1.3 gives O_epsilon(n^epsilon) completions for each coefficient triple (r1,r2,ell), since its kernel height is <=n. The exact input counts primitive integer solutions of AX^3+BY^3+CZ^3=0 of height <=H by O_epsilon(H^epsilon), uniformly in ALL nonzero signed integer coefficients. It does not count the frequency of coefficient triples.

Counting those triples by (6), and then recovering selected paths, gives
\[
 \#\{\text{paths}:\alpha\ge A,\gamma\ge C\}
 \ll_\epsilon n^\epsilon\frac{R_1^2R_2^2}{A^4C}.              \tag{8}
\]
If R1/A<1 there are no paths. Thus all floor/empty-range cases are covered. Starting with a smaller epsilon absorbs products of subpower factors.

Combining (8) with cycle continuation yields
\[
 \#\{C_6:\alpha\ge A,\gamma\ge C
             \text{ at some shared-right path}\}
 \ll_\epsilon n^\epsilon\frac{K_1R_1^2R_2^2}{A^4C}.           \tag{9}
\]
For Y>=1, partition alpha,gamma dyadically. On ranges [A,2A),[C,2C), alpha^4 gamma>=Y implies A^4 C>Y/32. At most O(log(n)^2) ranges occur. Hence
\[
 \#\{C_6:\alpha^4\gamma\ge Y
             \text{ at some shared-right path}\}
 \ll_\epsilon n^\epsilon\frac{K_1R_1^2R_2^2}{Y}.              \tag{10}
\]
The same argument, multiplying each dyadic contribution (9) by its maximal weight, proves
\[
 \sum_{C_6}\sum_{\text{shared-right paths}}\alpha^4\gamma
                 \ll_\epsilon n^\epsilon K_1R_1^2R_2^2.     \tag{11}
\]
Swapping summands gives these statements for shared-left paths, with R1,R2 interchanged and K1 replaced by K2. The two gcds always belong to the SAME path. No symmetric interchange of alpha,gamma on a fixed shared-right path is asserted.

## 4. The critical rectangle

Put
\[
 R_1=2n^{11/51}/(\kappa_1^2\kappa_2)^{1/3},\quad
 R_2=2n^{11/51}/(\kappa_1\kappa_2^2)^{1/3},\quad
 K_1=2n^{32/153}/(\kappa_1^4\kappa_2^2)^{1/9}.
\]
**Corollary 2 (Critical shaped sector).** The n^(26/51+epsilon) cycle budget holds in the sector (13), with the analogous coefficient-interchanged shared-left statement.

For sufficiently large n, (2) holds uniformly. The finitely many smaller nonempty targets are uniformly bounded since every positive parameter is bounded by a power of that fixed threshold. Exactly
\[
 K_1R_1^2R_2^2=
 32n^{164/153}/(\kappa_1^{22/9}\kappa_2^{20/9}).              \tag{12}
\]
Consequently the n^(26/51+epsilon) cycle budget is already met on paths with
\[
 \alpha^4\gamma\ge
 n^{86/153}/(\kappa_1^{22/9}\kappa_2^{20/9}).                 \tag{13}
\]
If this threshold is below one, (11) bounds all cycles since every path weight is >=1. For the shared-left orientation interchange the coefficient exponents. No squarefreeness of the coefficients is needed for this sector count.

For this specialization there is a second route: the existing [critical criterion, Section 4, Lemma 4](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md) supplies actual codegrees in precisely this shaped rectangle. Combining it directly with the coefficient/path count (8), the path-extension argument, and the dyadic sum proves the same critical sector. Thus the critical specialization does not require the separate general codegree argument in Section 2.

For coefficient one, the old individual alpha threshold is n^(43/306). The compatible exponent sizes alpha~n^(1/8), gamma~n^(2/25) satisfy
4/8+2/25=29/50>86/153,
1/8<43/306,
and 1/8+2/25=41/200<11/51.
This compares allowed sizes only; actual cycles in those ranges are NOT supplied.

Combining alpha*gamma<R2 with (13), every coefficient-one path in this sector must have alpha >> n^(53/459). This is necessary, not sufficient. In the balanced max box R1=R2=K1=T, (10) is O_epsilon(n^epsilon*T^5/Y), meeting the T^2 budget for Y>=T^3. Compatibility then requires alpha >> T^(2/3), compared with the old individual threshold T^(3/4).

## 5. Exact boundaries and an elementary obstruction

This is an actual-frequency bound: the new C^-1 factor comes from summing the large-square-divisor condition on ell BEFORE applying the coefficient-uniform primitive cubic theorem. It is not an assertion that O(1) completions at each slice bound the number of slices.

The proof does not show that every residual cycle has a large alpha^4 gamma path. It does not close the unrestricted critical-rectangle hypothesis or improve the proved global prescribed-sum exponent 3/10. Uniformity for each fixed coefficient pair does not authorize summing an unrestricted number of coefficient pairs. This is not a fixed-difference result or a global triple count.

The nonzero and squarefree hypotheses are essential in the two uses of prime-valuation parity. No zero-coefficient cubic is passed to Heath-Brown. The count uses one selected actual representation per edge, so it is a graph-cycle statement rather than an uncontrolled representation multiplicity statement.

The lifted factorization alone cannot strengthen alpha*gamma<R2 by a power. For odd coprime r>s, put
\[
 s_1=(r^2+s^2)/2,\quad s_2=(r^2-s^2)/2,\quad \alpha=rs.
\]
These are primitive integers with s1^2-s2^2=alpha^2 and alpha/s1=2rs/(r^2+s^2) tending to one along r/s->1. Multiplying s_i by a coprime gamma gives the same lifted pattern with alpha*gamma approaching the root cap. This is sharpness ONLY for the integer constraint, not an actual kernel-cube cycle or prescribed target.

## 6. Source and comparison scope

The [critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md), Section 4, supplies the comparison single-gcd bound and the optional shaped-codegree route. General codegree recovery and the joint square-divisor frequency argument needed for the standalone theorem are written above.

The sole non-elementary counting input is D. R. Heath-Brown, *Counting Square-full Solutions to x+y=z*, arXiv:2601.07817v1, Theorem 1.3:
https://arxiv.org/abs/2601.07817v1
The actual retained primary statement was reread. Its uniform primitive signed-cubic scope is stated exactly in Section 3. The global triple theorem in its introduction is not used.


## 7. Square-quotient strengthening: setup

Use precisely the actual selected graph (1), caps (2), path normalization (3), and nonzero integer lift (4) above. Thus the coefficients are fixed positive integers, roots and kernels are positive, kernels are squarefree, summands are coprime, and
\[
 1\le R_i,K_i\le n,\qquad 2(R_1R_2)^2<n.
\]
At a shared-right path with distinct left kernels put
\[
 e=|\ell_0|=\frac{|s_1^2-s_2^2|}{\alpha^2}
           =\frac{|c_1^2-c_2^2|}{\alpha^2\gamma^2},
 \qquad W=\frac{\alpha^3}{\sqrt e}.                         \tag{14}
\]
The square root is a positive REAL square root of a positive integer, used only as a counting weight, not as an element of a coefficient field. Here e>=1 is an integer. All counts are uniform in each fixed positive coefficient pair, and never an aggregate over varying pairs.

## 8. Counting small square quotients

**Theorem 3 (Small-quotient path frequency).** For real A,E>=1 the number of actual selected shared-right paths with alpha>=A and e<=E is
\[
 O_\epsilon\left(n^\epsilon R_1^2R_2\sqrt E/A^3\right).
                                                               \tag{15}
\]
If R1<A or R2<=A there are no such paths. Enlarging E beyond its actual range is harmless.

Indeed r1,r2<=R1/A and
\[
       0<|\ell|=\gamma^2 e\le L=R_2^2/A^2 .
\]
Count DISTINCT possible nonzero signed integers ell by the union bound
\[
 \#\{\ell\}\le
 2\sum_{1\le t\le E}\left\lfloor\sqrt{L/t}\right\rfloor
 \le 2\sqrt L\sum_{1\le t\le E}t^{-1/2}
 \le 4\sqrt{LE}.                                             \tag{16}
\]
The sums run over integers. If L<1 there is no ell; terms with t>L have zero floor. The estimate is valid for every real E>=1. Different square-quotient presentations of the same ell may be counted more than once, which only enlarges the upper bound. No claim of a unique presentation is needed.

There are therefore O(R1^2 R2 sqrt(E)/A^3) possible triples (r1,r2,ell). The exact cubic is still
\[
 \kappa_1r_1^2b_1^3-\kappa_1r_2^2b_2^3+\kappa_2\ell d^3=0.
\]
Its three coefficients are nonzero and its kernel triple is primitive because each b_i is coprime to d. The SAME coefficient-uniform Heath-Brown Theorem1.3 input gives n^epsilon completions per triple. A kernel triple determines at most one selected actual path. Thus (15) counts path FREQUENCY after counting the coefficient set; it does not substitute completion bounds for frequency. No value of alpha,gamma or e has to be added as an independent multiplicity after the selected kernel path is recovered.

The nonzero lift uses the actual squarefree distinct kernels: c1=c2 would imply a1^2 b1^3=a2^2 b2^3 and hence b1=b2 by prime-valuation parity. Coprimality gives gcd(alpha,gamma)=1 and makes ell0 integral. These are precisely the old gates, all retained.

## 9. Weighted path and cycle sums

**Theorem 4 (Square-quotient weighted cycle tail).** Under the same hypotheses,
\[
 \sum_{\text{selected shared-right paths}} W
                  \ll_\epsilon n^\epsilon R_1^2R_2,           \tag{17}
\]
and
\[
 \sum_{C_6}\ \sum_{\text{shared-right paths in }C_6}W
                  \ll_\epsilon n^\epsilon K_1R_1^2R_2.       \tag{18}
\]
Consequently for EVERY real Y>0,
\[
 \#\{C_6:\ W\ge Y\text{ at some shared-right path}\}
             \ll_\epsilon n^\epsilon K_1R_1^2R_2/Y.           \tag{19}
\]
There is a side-swapped version throughout.

For (17), partition alpha into dyadic ranges [A,2A) and e into [D,2D). Since 1<=alpha<=R1<=n and 1<=e<R2^2<=n^2, there are O(log(n)^2) pairs of ranges. Apply (15) with E=2D. The weight on that range is at most 8A^3/sqrt(D), so its total is O(n^epsilon R1^2R2), independently of A,D. Absorb the finite logarithmic factor by starting with a smaller epsilon.

For (18), a specified shared-right path has O(K1*n^epsilon) actual simple-cycle continuations. This is the unchanged two-codegree argument: choose a third left kernel, then its two common right neighbors. The general actual codegree proof uses the unit fraction (a1*c2)/(a2*c1), height<=R1R2, the bound2(R1R2)^2<n, and unique root/kernel recovery; it was proved in Section2. Multiplying (17) by the continuation bound proves (18).

Every counted cycle in (19) contributes at least Y to (18). Positivity of W proves (19) directly. W NEED NOT be an integer or be at least one. Its tail parameter is Y>0, not only Y>=1.

## 10. The original joint-gcd estimate is a corollary

The lifted equation gives the STRICT pointwise inequality
\[
 \alpha\gamma\sqrt e
       =\gamma\sqrt{|s_1^2-s_2^2|}<R_2,
 \qquad \alpha^4\gamma<R_2W.                                \tag{20}
\]
Thus the weighted bound and tail in the original Theorem1 follow immediately from (18)--(19), with prefactor K1 R1^2 R2^2. A direct independent old route also remains valid: count ell with a square divisor gamma^2>=C^2 as in the original Sections1--3.

The old two-threshold path bound is recovered as well. If alpha>=A and gamma>=C, then e<R2^2/(A^2 C^2). If this upper bound is at most one the paths are empty. Otherwise (15) gives
\[
 \#\{\text{paths}:\alpha\ge A,\gamma\ge C\}
        \ll_\epsilon n^\epsilon R_1^2R_2^2/(A^4C).
\]
Hence retaining the earlier result as a corollary loses no prior statement.

The new sector is larger at the level of allowable sizes. For coefficient one, alpha~n^(1/8), gamma~1 and e~1 satisfy W~n^(3/8)>n^(53/153), whereas alpha^4 gamma~n^(1/2)<n^(86/153). The lifted cap alpha*gamma*sqrt(e)<R2 is compatible with these exponents. This does NOT construct actual kernel cubes, actual paths or cycles of those sizes.

## 11. Critical and balanced specializations

For the shaped critical caps,
\[
 K_1R_1^2R_2
  =16n^{131/153}/(\kappa_1^{19/9}\kappa_2^{14/9}).             \tag{21}
\]
Therefore the n^(26/51+epsilon) cycle budget is proved on the sector
\[
 W\ge
 n^{53/153}/(\kappa_1^{19/9}\kappa_2^{14/9}),                 \tag{22}
\]
with coefficient exponents interchanged for shared-left paths. The coefficient-one sector necessarily has alpha>=n^(53/459)*e^(1/6), by rearranging W=alpha^3/sqrt(e)>=n^(53/153).

There is an alternative route for this SPECIFIC specialization: apply (15), then use the existing shaped actual-codegree fact from the critical criterion Section4 Lemma4 instead of the local general codegree theorem, and repeat the dyadic weighted-sum proof. This never broadens that existing fact to arbitrary rectangles.

For all actual paths,
\[
                         W>1/R_2,                          \tag{23}
\]
by (20) and alpha^4 gamma>=1. Thus (22) controls ALL cycles if its threshold is <=1/R2. A threshold merely below one is insufficient for that conclusion. The sector statement itself is valid for every positive threshold by (19), including thresholds between1/R2 and1.

In the balanced maximum box R1=R2=K1=T, the new tail is O(n^epsilon*T^4/Y). The T^2*n^epsilon budget therefore holds on W>=T^2. This contains the previous alpha^4 gamma>=T^3 sector by (20), and permits the e-small ranges not selected by that old condition.

## 12. Source overlap, scope and remaining obligation

No new imported theorem is used. The sole non-elementary input remains the EXACT Theorem1.3 of Heath-Brown2601.07817v1 for primitive integer diagonal-cubic zeros, uniformly in NONZERO SIGNED coefficients. The ordinary codegree proof and selected-path recovery are unchanged. Equation(16), not the completion theorem alone, supplies the additional coefficient-frequency restriction.

The small-square-quotient count(16) refines the coefficient frequency used in(8), with the same uniform cubic input.

The stronger tail still does not force every residual cycle to have a large W path. It does not prove the full critical rectangular C6 hypothesis, its residual equivalent, a global prescribed-sum exponent below the proved3/10, or a uniform varying-coefficient aggregate. Small W actual frequency remains unbounded by this argument.
