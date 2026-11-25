# Dual weighted-difference codegrees for varying shifts

## A finite graph bound

Fix \(C,\varepsilon>0\). Let \(X\ge2\), fix positive integers \(F,G\), and fix a sign \(\eta\). Consider any collection of distinct triples of labels (x,y,k) in positive integers with

    T = F x y + eta G k > 0,            gcd(y,k)=1.

Assume F,G,x,y,k,T are at most X^C. The x-values range over a set of cardinality A and the possible column pairs (y,k) over a dictionary of cardinality B. For each selected label triple choose one representation

    T = W z^q,       W,z positive integers,       q >= 2 integer.

One can require q <= max(2,C log_2 X) without losing any representation when the exponent is free, by using q=2 if z=1. This normalization increases a prescribed (W,q)-dictionary by at most a factor of two and preserves a W-only dictionary. Alternatively retain any prescribed finite set of exponent choices fixed before the height limit. Suppose A_W is the number of distinct row labels \((x,W,q)\) after the choice, and M is the number of edges. Then

    M <<_{C,epsilon} X^epsilon (A_W + B sqrt(A_W)).                  (D1)

In particular A_W <= A times the number of records (W,q). When the number of q choices is logarithmic or bounded, one may use the number of W choices in this product, with the resulting loss absorbed in X^epsilon. The normalized coefficient and record caps are bounds on the actual dictionary. No hypothesis of independence of x and W is used.

Use a bipartite graph with row label (x,W,q) and column label (y,k), selecting the edge from the chosen representation. It is simple since the labels determine x,y,k, and every original label triple is counted once. No complete-value uniqueness across different triples is required. Two distinct columns (y_1,k_1),(y_2,k_2) have distinct ratios k_j/y_j, since each positive pair is primitive. At a common row, multiply the first edge equation by y_2 and the second by y_1, then subtract:

    W (y_2 z_1^q - y_1 z_2^q) = eta G (y_2 k_1 - y_1 k_2) = D.

The integer D is nonzero, has height at most 2X^(3C), and is fixed by the two columns. Thus W divides |D|, with divisor-many choices independent of the target coefficient dictionary. For W and q fixed, the right side D/W is fixed and nonzero. A solution (z_1,z_2) fixes x from either edge equation. For any fixed bound Q on q, Zeng's Lemma 3.1 gives O_{Q,C,epsilon}(X^epsilon) such root pairs, uniformly in the two positive weights and in D/W. This lemma requires only polynomial heights and no coprimality between the two roots or the weights. Roots here satisfy z_j <= X^(C/q). For q>Q, choosing z_1 gives at most X^(C/Q) choices, and then z_2 is unique if it exists. Set a fixed Q>10C/epsilon and sum at most O_C(log X) exponent choices on the normalized labels, along with the divisor count on W. Taking the small epsilon losses at the individual stages proves that the column-pair codegree is O_{C,epsilon}(X^epsilon), uniformly even with exponent choices varying with X within the normalized range. Arbitrary finite choices unnormalized but fixed in X can instead be summed with a constant depending on that set. The reduced positive pairs exclude \(D=0\). Without that hypothesis, coincident slopes can have large codegree.

Summing the numbers of pairs of neighbor columns over rows gives at most O(X^epsilon B^2). Cauchy–Schwarz then gives M^2/A_W - M << X^epsilon B^2. Solving this inequality and reducing auxiliary epsilon proves (D1). Its linear term matters if the selected graph is too sparse to give repeated edges in a row. If original variable labels have redundant representations, choose one representation per counted object before this graph argument or bound the multiplicity separately. In an anatomic box at fixed depth, divisor-many monomial representations can be absorbed in X^epsilon.

In the [varying-shift notation](ABC_AGGREGATE_SHIFT_DICTIONARY_STARS.md), put n=X^L. Suppose the complete unfrozen U factor has root cap a, x uses cap p, y uses cap a-p, the shift dictionary k has cap b, and the target coefficient record (W,q) has cap f. At each fixed L, allow arbitrary strictly smaller margins to absorb box losses. Then A_W <= n^(p+f+o(1)) and B <= n^(a-p+b+o(1)). Consequently the fibre exponent r satisfies

    r <= max{p+f, a+b+f/2-p/2}.                                  (D2)

No numerical height gap is needed and the shift dictionary is allowed to overlap the whole numerical range. If the given decomposition admits any chosen split p and p=(2(a+b)-f)/3 is available in [0,a], balancing the two terms gives r <= 2(a+b+f)/3. This last optimization presupposes the specified row block exists. The inequality itself holds for each available p. Fixed restrictions on the shift signs or coefficient height only enter at the initial graph. The result complements the existing modulus-weight uniqueness criteria, whose height assumptions are different. It adds no unconditioned estimate on the dictionary cap f of a fully diffuse target.

## Complete-value graph at the minimum nonatomic endpoint

Let A,B be dictionaries of complete values and count a set of distinct additive edges A+B=T with each target supplied one power record W z^q, q>=2. Use x=A, y=1, k=B and F=G=1 above. Then D=B_1-B_2 and each weighted difference is unweighted after dividing by W. There is also a direct coefficient-independent proof here: for positive roots v,w, if v^q-w^q=m!=0, the signed difference v-w is a nonzero divisor of m. At a fixed signed difference, the left side is strictly monotone in the smaller root, so there is at most one pair. This argument is uniform in q>=2 and summing q still costs at most a logarithm after normalizing q for root one. There is no reliance on a new number-field estimate in this special case. The bound is invariant under changing the sign of B or interchanging A and B.

Use the existing complete-label subgraph at d=1/3, total nonatomic mass Gamma=1/2 and c=1/6, with n=X^L. The two free complete dictionaries have capacities c+o(1) and the edges approach exponent 2c. This statement is used on the extracted complete-label subgraphs, after taking the selection losses sufficiently small relative to the fixed margin in the next inequality. A target dictionary of coefficient and exponent records with normalized exponent u yields

    count exponent <= max{c+u, 3c/2+u/2}.

At the endpoint one must have u >= c=1/6. More quantitatively, for any delta in (0,c), the edges using an arbitrary subdictionary of size at most n^(c-delta) number at most n^(2c-delta/2+o(1)). Relative to a subgraph of edge exponent 2c-o(1), this is an exponentially small fraction in log n. The same gives a lower tail bound for the information of the selected power record, since at most n^(c-delta) records can have edge probability at least n^(-c+delta). If an anatomic decomposition collects all layers divisible by q into the q-th-power root, the complementary radical mass bounds u, because the coefficient is determined by the remaining layers and the fixed anatomy coefficient. Thus at a complete-label minimum endpoint the q-divisible mass on each letter cannot exceed c by a fixed margin. This statement allows q to vary with the outer small budget, and uses the height limit first.

This sharpens the [dictionary threshold from common powers](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-value-codegrees-from-common-powers), where (FD1) charges two freely varying coefficients and only yields u>=1/12 for singleton complete labels. It does not duplicate the mixed-power row codegrees (MP1) or (MS8), which fix several complete rows and count the common columns by higher-degree curves. The new bound classifies the rows by their chosen record and applies a difference to two fixed columns sharing such a row.

## An exact numerical separation of the dictionary inequalities

Set c=1/6 and use a target cofactor cap u=1/8. Then (FD1) with singleton complete labels allows it, since 2u=1/4>c. The target exponent caps in (MS7) and (MS8) also allow it whenever a large power is admissible, because m u>=3/8 for m>=3. In contrast (D1) gives

    max(c+u,3c/2+u/2) = max(7/24,5/16) = 5/16 = 1/3-1/48.

This cofactor cap is compatible with the following formal valuation data. Take an integer K tending to infinity, L=K^(-4) and q=K^3. On the target put mass 1/6 at original valuation 2K^4. Distribute total mass 1/24 equally over the values qj, j=3K,...,5K, and mass 1/8 equally over qj-1 and qj+1 at the same j. The total mass is 1/3, the first moment in scaled degrees is 1, and the limiting nonatomic mass is 1/6 distributed uniformly over [3,5]. The first two parts are q-divisible and the remaining part is not, so its complementary mass is precisely 1/8. Here q is fixed before each height limit, Lq tends to zero and Lq^(3/2) tends to infinity. The example is a numerical comparison for these specific dictionary inequalities and does not assert viability of the profile under the entire arithmetic program.

## Cofactor ratio dictionaries for two fixed columns

There is a second way of using the same nonzero determinant D, which does not add W to the row label. Start with the graph on x and (y,k) and one power representation chosen per edge. For each pair j,j' of distinct columns, and fixed exponent q, let R_{j,j'}(q) be the number of W_j/W_j' ratios occurring over those common rows for which both edges have exponent q. Write a ratio in lowest terms as r/s. Then W_j=g r and W_j'=g s for a positive integer g. Eliminating x gives

    g (y_{j'} r z_j^q - y_j s z_{j'}^q) = D.

Thus g divides |D| and every fixed reduced ratio gives only X^epsilon common rows, by the same weighted-difference lemma and the same elementary treatment of large q. The root pair determines x. All integer factors in this equation have polynomial height, since the original records did. In particular, no independence of cofactor values on the two edges is needed. Positivity and column primitivity exclude D=0 just as before.

With a single fixed q, if M is the edge count, then

    M^2/A - M <<_{C,epsilon} X^epsilon sum_{j!=j'} R_{j,j'}(q).       (D3)

For m possible exponents, use the graph of each exponent separately and sum. By Cauchy–Schwarz over the m parts, the left side becomes M^2/(m A)-M and the right side gains a sum over q. Normalization bounds m by O(log X), so this has the same exponent implications. Alternatively fix one exponent on a largest part before using (D3).

At complete-label capacities n^(c+o(1)) on both free letters with edge count n^(2c-o(1)), there must be n^(2c-o(1)) column pairs with ratio dictionary of size n^(c-o(1)), after retaining one exponent. More explicitly fix epsilon_0>0 and take the capacity and count losses delta sufficiently small relative to epsilon_0. The sum of common-row counts over column pairs has exponent at least 3c-3delta, while pairs with fewer than n^(c-epsilon_0) common rows give at most exponent 3c+2delta-epsilon_0. There remain at least n^(2c-4delta+o(1)) pairs, because each codegree is at most n^(c+delta). Every retained pair then gives at least n^(c-epsilon_0-o(1)) distinct cofactor ratios by the fixed-ratio bound. Let the losses decrease in the endpoint sense. This is stronger than a lower bound on the cardinality of the whole coefficient dictionary when the ratios have special structure.

For instance, if W_{ij}=U_i V_j on the selected edges for positive rational row and column factors, then each column pair has at most one cofactor ratio, regardless of the cardinality of the coefficient dictionary. Equation (D3) gives M << X^epsilon(A+B sqrt(A)). This would save up to c/2 at the complete-label endpoint. More generally any uniform ratio cap V gives M << X^epsilon(A+B sqrt(A V)). This concerns the ratios of the target power cofactors, which are distinct from the complete A/C cross ratios in (XR1).

## Sources

- Chunping Zeng, *Complete modulus locking in the exceptional set of the abc conjecture*, Hexagon `2611.00621v1`, §3.1, Lemma 3.1, https://hexagonmath.org/pdf/2611.00621v1. The analytic lemma is used independently of the numerical certificate.
- Hongmin Li, *Squares and twisted energy in the exceptional set of the abc conjecture*, rStation math paper 28, version 4, §3.1 Remark 3.2a and Appendix B, https://rstation.ai/problems/math/papers/28/v4. This is the antecedent arbitrary-fixed-degree weighted-difference bound reproduced by Zeng.
- Ernie Croot and Chi Hoi Yip, *Diophantine tuples and product sets in shifted powers*, JLMS 113 (2026), e70499, Proposition 4.4, https://doi.org/10.1112/jlms.70499. This is a close graph analogue for a single fixed shift with both factors at least 2|shift|^17. That proposition is not needed for (D1).
- [Complete values determined by single outside labels](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-values-determined-by-single-outside-labels) and [reduction close to minimum mass](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass).
