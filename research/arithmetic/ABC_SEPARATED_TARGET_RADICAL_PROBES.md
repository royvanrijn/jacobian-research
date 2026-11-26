# Private-prime packing for multiplicatively separated target probes

## Setup and pointwise proof input

Let \(X/2<c\le X\), \(a+b=c\), with positive primitive integers. Put \(n=X^L\), where \(L>0\) is sufficiently small and fixed before \(X\) tends to infinity. Fix two different letter roles, denoted \(U,V\), and let \(W\) be the remaining letter. On a fixed choice of roles one has \(W=\epsilon U+\epsilon' V\) with fixed \(\epsilon,\epsilon'\in\{-1,1\}\). Let \(\mathcal E\) be any graph of actual pairs \((U,V)\), with \(H_U,H_V\) complete labels and no duplicate edges. Keep all three raw radical caps
\[
\operatorname{rad} T\le n^{r_T},\qquad T\in\{U,V,W\},\qquad 0\le r_T\le C
\]
and any other conditions defining the original cell. The conclusions also apply on a subgraph of a residue cell, with all the original caps. Write \(R(U,V)=\operatorname{rad}(W(U,V))\).

Use the fixed-letter bound FL1. [Fixed-letter radical caps](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools) For any fixed \(B,C,\eta>0\), sufficiently small fixed \(L>0\) and any \(\varepsilon>0\), it bounds a fibre of one fixed actual letter by
\[
O_{B,C,\eta,L,\varepsilon}\bigl(n^{\min(\ell_1,\ell_2)/2+\eta}X^\varepsilon\bigr)
\]
when the two free letters have outside radical exponents \(\ell_1,\ell_2\in[0,C]\) after prescribed prime pools of integer product at most \(X^B\). In our use one free cap is its original cap, with empty pool. The other is the minimum of its original cap and the derived outside cap. All constants in the estimates below are uniform in original letters, residues, model functions and the actual pools. A derived outside exponent of zero can also use the elementary fixed-support bound directly.

## SRP1. Rectangle denominator packing

For a nondegenerate rectangle with anchors \(U_1,U_2,V_0\) and final query \(V\), write \(R_{i0}=R(U_i,V_0)\), \(R_{iV}=R(U_i,V)\) and
\[
D(V)=\frac{R_{10}R_{2V}}{R_{1V}R_{20}}.
\]
At a prime \(p\mid W(U_1,V)\) off the fixed pool of prime factors of \(W(U_1,V_0)(U_1-U_2)\), one has \(p\nmid R_{10}\) and \(p\nmid R_{2V}\), since \(W(U_1,V)-W(U_2,V)=\epsilon(U_1-U_2)\ne0\). The prime \(p\) therefore divides the reduced denominator of \(D(V)\). Consequently, a bound \(n^\beta\) on the radical of that denominator forces the outside radical of \(W(U_1,V)\) to be at most \(n^\beta\). The pool has integer product at most \(X^2\). Apply (FL1) with fixed \(U_1\). The number of such \(V\) for each realized anchor is
\[
O\bigl(n^{\beta/2+\eta}X^\varepsilon\bigr),
\]
and the number of such rectangles is \(O(H_U^2H_V n^{\beta/2+\eta}X^\varepsilon)\). The estimate applies a fortiori with bounded rational height of \(D(V)\), and reversing \(U_1,U_2\) gives the numerator version. A reference query \(V_0\) and a candidate \(V=V_0\) can be included at cost one.

If \(H_U,H_V\le n^{u+\tau}\) and \(|\mathcal E|\ge n^{2u-\tau}\) with \(u>0\) fixed and \(0<\tau\ll u\), there are at least \(n^{4u-O(\tau)}\) nondegenerate ordered rectangles for large \(X\). Indeed the rectangle homomorphism count is at least \(|\mathcal E|^4/(H_U^2H_V^2)\) by two applications of Cauchy-Schwarz, and repeats cost at most \(|\mathcal E|(H_U+H_V)\). Thus rectangles with reduced denominator radical at most \(n^{2u-\xi}\) are power-small whenever \(\xi>0\) is fixed and \(\tau,\eta,\varepsilon/L\) are sufficiently small in terms of \(u,\xi\).

## SRP2. An arbitrary multiplicative model of bounded rank

Fix \(r\ge0\) and suppose that on every actual edge of a subgraph there is an exact rational identity
\[
R(U,V)=C(U,V) A(U)B(V)
               \prod_{j=1}^r F_j(V)^{k_j(U)}.                     \tag{1}
\]
The one-label functions \(A,B,F_j\) are arbitrary positive rationals. Exponents \(k_j(U)\) are arbitrary integers. Their sizes and the heights of the one-label functions are unrestricted. Assume \(C(U,V)\) is a positive divisor of \(R(U,V)\) and \(C(U,V)\le n^\beta\). Let \(m=r+2\). One may instead use \(m=d+1\) if the matrix with rows \(1,k_1,\dots,k_r\) has rank at most \(d\) on every anchor list. If \(B\) is identically 1, the row of ones can be omitted. The case \(m=1\) after omitting that row is the direct one-anchor fixed-fibre bound below, with any realized \(V_0\) or no query to count.

For any \(m\) distinct \(U_i\), the matrix of \(m\) columns has a nonzero integral null vector \((\lambda_i)\). Choose a realized common \(V_0\) if one exists. Divide (1) at \((U_i,V)\) by (1) at \((U_i,V_0)\), and take the product with respective powers \(\lambda_i\). Both \(A(U_i)\) and the factors in \(V,V_0\) disappear. Therefore
\[
\prod_i\left(\frac{R(U_i,V)}{R(U_i,V_0)}\right)^{\lambda_i}
 =
\prod_i\left(\frac{C(U_i,V)}{C(U_i,V_0)}\right)^{\lambda_i}.        \tag{2}
\]
Freeze an \(i\) with \(\lambda_i\ne0\). Use the fixed pool from the nonzero differences \(U_j-U_h\) and all \(W(U_j,V_0)\). Its prime product is at most \(X^{\binom m2+m}\). A prime \(p\mid W(U_i,V)\) outside this pool occurs in exactly one of the \(W(U_j,V)\), and no reference letter. The left side of (2) then has \(p\)-adic valuation \(\lambda_i\). The containment \(C(U_j,V)\mid R(U_j,V)\) forces \(p\mid C(U_i,V)\), because otherwise the right side valuation vanishes. In fact the positive valuation of that residual is 1. The outside radical of the actual letter \(W(U_i,V)\) is thus at most \(n^\beta\). With \(U_i\) fixed, (FL1) bounds the entire common neighbourhood by
\[
K_m=O_{r,C,\eta,L,\varepsilon}\bigl(n^{\beta/2+\eta}X^\varepsilon\bigr).
\]
Nothing is paid for the reference letter, which is merely chosen in the proof of a uniform bound for each anchor list. No enumeration of exponent vectors or model functions is used. If \(C\equiv1\), predetermined fixed-support counting gives \(K_m\ll_{\varepsilon,m}X^\varepsilon\).

Write \(d(V)\) for the degree of \(V\). There are at most \(K_m H_U^m\) ordered \(m\)-stars by summing over ordered distinct anchor lists. The vertices with \(d(V)<2m\) give at most \(2mH_V\) edges. For the others, \((d(V))_m\ge(d(V)/2)^m\), and the power mean inequality gives
\[
|\mathcal E|\ll_m H_V+
       H_U H_V^{1-1/m}K_m^{1/m}.                                  \tag{3}
\]
In particular if \(H_U,H_V\le n^{u+\tau}\), \(u>0\), the exponent is bounded by
\[
(2-1/m)u+\frac{\beta}{2m}+(2-1/m)\tau+\eta
\]
with a factor \(X^\varepsilon\), taking inputs at suitably reduced \(\varepsilon,\eta\). For fixed \(r\) this gives a strict saving for \(\beta<2u\) with a fixed margin and sufficiently small \(\tau,\eta,\varepsilon/L\). If there are at most \(n^{b+o(1)}\) model charts with those common bounds, their catalogue costs \(b\) in the final exponent. The saving remains if \(mb+\beta/2<u\) by a fixed margin and the errors are small. For exact separation \(R=A(U)B(V)\), use \(m=2\) and obtain \(|\mathcal E|\ll n^{3u/2+O(\tau)}X^\varepsilon\). The \(m=2\) argument also proves the approximate-separation exponent \(3u/2+\beta/4\) up to losses. Any separation of the full triple radical with the same shape gives a separation of \(R\), by division by the one-letter factors \(\operatorname{rad}(U)\operatorname{rad}(V)\).

The same exact bound holds if each \(C(U,V)\) is a rational with prime support contained in \(W(U,V)\) and its support product is at most \(n^\beta\). Without support containment, the proof instead requires that, for any anchor list, the union of the supports of the active residual factors \(C(U_j,V)\) with \(\lambda_j\ne0\) has product at most \(n^\beta\), and adds the reference residual supports to the fixed pool. Polynomial height of that augmented pool suffices, and follows from fixed \(m,\beta\) and the stated support bound on actual queries. A cap \(n^{\beta_0}\) on each residual support then gives \(\beta=m\beta_0\).

## SRP3. Modular depth and simultaneous probes

Suppose the factors on the right of (1) and \(R\) are units in a common fixed \(\mathbb Z/M\mathbb Z\), and (1) is initially a congruence. Keep \(C(U,V)\mid R(U,V)\). For each anchor list assume the existence of a nonzero integral null vector with
\[
M> n^{r_W\|\lambda\|_1}.                                           \tag{4}
\]
The congruence product, applying \(\lambda_i\), now says
\[
\prod_i\left(\frac{(R/C)(U_i,V)}{(R/C)(U_i,V_0)}\right)^{\lambda_i}
        \equiv1\pmod M.
\]
Cross-multiplication uses only the actual positive integers \(R/C\), each at most \(n^{r_W}\). Both cross-products are at most \(n^{r_W\|\lambda\|_1}\), so (4) makes their equality exact, giving (2). The earlier proof then applies. Neither the one-label functions' numerical heights nor their representatives in the residue ring enter (4).

For example if \(|k_j(U)|\le K\) and \(K\ge1\), the null vector can be chosen with
\(\|\lambda\|_1\le (r+2)(r+1)!K^{r+1}\). For this, choose independent columns and equally many independent rows giving a nonsingular minor of the matrix, add one extra column and take the resulting signed maximal minors. The relation on these rows annihilates all rows by their rank. Fixed \(r,K\) and \(M\ge X^\nu\) with fixed \(\nu>0\) imply (4) for sufficiently small fixed \(L\). With \(K\) varying with \(L\), one needs the explicit comparison \(Lr_W(r+2)(r+1)!K^{r+1}<\nu\). The exact model required no such bound. In particular, exact separation modulo a depth-\(X^\nu\) modulus lifts at the rectangle level under \(2r_WL<\nu\). There is no assertion that an arbitrary congruence model with unbounded weights lifts.

## Simultaneous target probes and lost support

The same argument treats finitely many positive integer probes \(P_h(U,V)\) supported on the primes of their actual \(W(U,V)\). Suppose that for each \(h\) there is an exact identity
\[
P_h(U,V)=C_h(U,V) A_h(U)B_h(V)
             \prod_{j=1}^{r_h}F_{hj}(V)^{k_{hj}(U)}.
\]
Assume each residual \(C_h(U,V)\) is a divisor of its probe. On an actual edge, require that the union of the residual supports has product at most \(n^\beta\), and that the primes dividing \(W(U,V)\) missed by all probes have product at most \(n^\kappa\). Choose a single null vector for the row of ones and all weight rows from all probes, using \(m=2+\sum_h r_h\). The omitted row of ones can be dropped only when all the corresponding \(B_h\) are identically 1. For the one chosen active \(i\), any outside prime seen by a probe has left-hand valuation \(\lambda_i v_p(P_h(U_i,V))\ne0\) in that probe's identity. It therefore occurs in its corresponding residual. The outside radical for actual \(W(U_i,V)\) is at most \(n^{\beta+\kappa}\), giving (3) with this new price. These statements do not require squarefree probes. For congruence models in each probe, the depth comparison uses a true numerical height bound \(P_h/C_h\le X^{h_h}\) and requires \(M_h>X^{h_h\|\lambda\|_1}\) for each corresponding modulus. A radical bound on \(P_h\) alone does not provide that numerical height bound.

## SRP4. Exact models in both orientations

There is also a stronger two-anchor bound under an exact identity in which the weights are on \(V\), specifically
\[
R(U,V)=C(U,V) A(U)B(V)\prod_{j=1}^r G_j(U)^{k_j(V)}.
\]
Fix distinct \(U_1,U_2\). The valuation vectors of the realized positive rational ratios
\[
Q(V)=\frac{(R/C)(U_1,V)}{(R/C)(U_2,V)}
\]
lie in one affine rational space of dimension at most \(r\), since the ratio of the edge models cancels \(B(V)\) and uses the fixed generators \(G_j(U_1)/G_j(U_2)\). Choose at most \(r+1\) realized ratio vectors spanning their affine hull. Any prime missed by all these ratios is missed by every ratio, as that prime coordinate vanishes on their affine hull. The product of supports of all selected ratios is at most \(X^{2(r+1)}\) since they are ratios of actual \(R/C\)'s, each at most \(X\). Adjoin the prime factors of \(U_1-U_2\). If \(p\mid W(U_1,V)\) is outside this pool, it does not divide \(W(U_2,V)\), and the zero \(p\)-valuation of \(Q(V)\) forces \(p\mid C(U_1,V)\). This gives \(K_2\ll n^{\beta/2+\eta}X^\varepsilon\) and the stronger \(m=2\) version of (3) for each fixed \(r\). It uses only the heights of realized post-residual ratios, so the generators need no height bound, and the integer weights remain unbounded. This affine-support argument concerns exact identities, and does not itself convert a congruence model. The original \(m=r+2\) method with \(U,V\) interchanged does extend to this orientation when the modular depth comparison holds.


Both orientations can also be present together. If the
model for \(R/C\) is
\[
 A(U)B(V)
    \prod_{\ell=1}^r F_\ell(V)^{k_\ell(U)}
    \prod_{j=1}^s G_j(U)^{h_j(V)},
\]
fix \(m=r+2\) anchors and a nonzero integral vector
\(\lambda\) killing the row of ones and all the anchor
weight rows \(k_\ell\). The valuation vectors of
\(\prod_i (R/C)(U_i,V)^{\lambda_i}\) lie in an affine
space of dimension at most \(s\) on the common queries.
Use at most \(s+1\) realized queries spanning that
space. Their products have support contained in the
union of the finitely many *actual* \(W(U_i,V_j)\)
supports, a pool of polynomial-height prime product
independent of the size of the weights. Adjoin all
pairwise anchor differences. Outside the pool the
valuation vector of this product is zero, and a prime
on an active \(W(U_i,V)\) can occur in no other anchor
role at that query. Thus it belongs to \(C(U_i,V)\).
The same fixed-anchor bound and (3) apply. Exchanging
orientations permits \(m=2+\min(r,s)\). The case of one
orientation only allows \(m=2\) by this proof and keeps
unbounded integer weights in exact identities.

A low-cost separate-label model is therefore excluded
on a full positive-fibre balanced cell if the target
radical left over from the specified probes has price
strictly less than \(2u\), subject to the stated model
catalogue costs. The argument is a fixed-original-fibre
bound and the graph retains all three original caps.
The question is whether actual cells force such models
or can be excluded by an estimate that permits the
nonseparated radical growth.

## SRP5. Fixed lock support and generic cross-ratio height

An external prescribed prime pool of polynomial prime product can be included in each proof above. In particular, replace \(R\) by the retained radical outside this pool, require \(C\) to divide the retained radical, and measure missed probe support only outside the pool. Such primes are adjoined to the pool in every use of (FL1), so their removal costs no residual price. In a congruence model arrange that the external pool contains the support of the common fixed modulus when needed. Then the retained \(R\) is a unit there, even when the modulus divides the actual target. The numerical depth comparison remains an upper bound on the actual retained probe, and all original raw caps stay in force.

There is a stronger generic cross-ratio conclusion, using the same private-support principle as PV1. [Private prime supports on actual boxes](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md) Suppose \(H_U,H_V\le n^{u+\tau}\) and \(|\mathcal E|\ge n^{2u-\tau}\). Fix \(u,\xi>0\) with \(\xi<2u\), fix integers \(m,k\ge2\), and take \(\tau,\eta,\varepsilon/L\) sufficiently small in terms of these parameters. In all but a power-small proportion of the realized injective ordered \(K_{m,k}\) arrays, every target position has private radical greater than \(n^{2u-\xi}\) outside the external pool. Privacy here means that the prime does not occur at any other complete target position and does not divide any complete row or column label.

For a distinguished target \(W(U_i,V_j)\), hold all \(U\) anchors and all but \(V_j\) fixed. Add to the pool the supports of all fixed target letters, all fixed complete row and column labels, and the nonzero differences between row anchors. Any prime of this target outside the pool is absent from every target in its column except this one, because sharing with another such target would imply divisibility of the corresponding row difference. It is absent from the moving column label \(V_j\) by primitivity of this actual edge, and is absent from all other complete labels by construction. If its private radical is at most \(n^{2u-\xi}\), the same is true for the outside radical in the augmented pool. By (FL1), fixing \(U_i\) gives \(O(n^{u-\xi/2+\eta}X^\varepsilon)\) possibilities for \(V_j\). Summing over the other anchors yields at most
\[
O\bigl(n^{(m+k-1)(u+\tau)+u-\xi/2+\eta}X^\varepsilon\bigr)
\]
arrays failing this test. By two applications of the power mean, the number of ordered array homomorphisms is at least \(H_U^m H_V^k p^{mk}\), where \(p=|\mathcal E|/(H_U H_V)\). The lower count and upper dictionary sizes give \(H_U,H_V\ge n^{u-2\tau}\) and \(p\ge n^{-3\tau}\), so the lower bound is \(n^{(m+k)u-O_{m,k}(\tau)}\). Repeated anchors cost at most \(O_{m,k}(n^{(m+k-1)(u+\tau)})\). The upper bound on failures is power-smaller after summing the \(mk\) tests.

On such a good realized rectangle the two reduced cross-ratio numerator factors have disjoint private radical portions, and neither portion appears in the denominator. The same holds for the denominator. Thus the numerator and denominator of \(D(V)\) each have radical greater than \(n^{4u-2\xi}\), improving the one-query denominator packing threshold for almost all rectangles. More generally, take any nonzero integer matrix \((b_{ij})\) and form the reduced rational \(\prod_{i,j}R(U_i,V_j)^{b_{ij}}\) in a good array. Its numerator is at least \(n^{(2u-\xi)\sum_{b_{ij}>0} b_{ij}}\), and the analogous denominator lower bound has \(\sum_{b_{ij}<0}|b_{ij}|\). These inequalities follow at the disjoint private primes, so hold simultaneously for all integer matrices, with no enumeration or bound on the weights. All bounds are on the actual prime supports of the original letters.

## Polynomially growing modular weights

The determinant bound on \(\lambda\) is unnecessary if more anchors are available. For a matrix with \(m>r+1\) columns, a row of ones and \(r\) signed integer weight rows bounded by \(K\ge1\), there is a nonzero integral null vector satisfying
\[
 \|\lambda\|_1\ll_{m,r} K^{r/(m-r-1)}.                             \tag{6}
\]
Indeed map \(\{0,\ldots,Q\}^m\) into its sum coordinate and \(r\) weighted sums. For \(Q\ge1\) the image has at most \((mQ+1)(2mKQ+1)^r\le(3m)^{r+1}K^r Q^{r+1}\) points. For \(Q\) a sufficiently large \(m,r\)-dependent multiple of \(K^{r/(m-r-1)}\), this is less than \((Q+1)^m\). Subtract a collision. Thus in the one-orientation modular result, if \(K(L)\le C_1 L^{-B}\) for fixed \(B,C_1\), choose fixed \(m>r+1+Br\). Then \(L\|\lambda\|_1\to0\), and any fixed positive numerical depth suffices. The exponent saving in (3) is still fixed when \(u-\beta/2\) has a fixed positive lower margin.

For both orientations of (5), suppose their integer weights satisfy \(K_U(L)\le C_1 L^{-B_U}\) and \(K_V(L)\le C_1 L^{-B_V}\), with \(r,s,B_U,B_V\) fixed. Let the common modulus satisfy \(M\ge X^\nu\) for fixed \(\nu>0\), and assume all units as before. Choose fixed integers \(m>r+1\), \(k>s+1\) with
\[
 \frac{B_U r}{m-r-1}+\frac{B_V s}{k-s-1}<1.
\]
On any injective actual \(K_{m,k}\), apply (6) independently to the row and column weight matrices. Their nonzero null vectors have \(L\|\lambda\|_1\|\mu\|_1\to0\) uniformly. At sufficiently small fixed \(L\), lift the tensor congruence exactly as before. At some position its coefficient \(\lambda_i\mu_j\) is nonzero. At a prime private to that position's target, the valuation identity forces that prime into the chosen residual. Thus every such array has at least one position with private radical at most \(n^\beta\), regardless of which positions are active in the two null vectors. This avoids imposing a nonzero coefficient at a prescribed query. The argument with multiple probes needs the same common null relations for all probes and the analogous depth checks, and replaces \(\beta\) with \(\beta+\kappa\).

Consequently, for \(H_U,H_V\le n^{u+\tau}\), the number of injective ordered \(K_{m,k}\) arrays is
\[
O\bigl(n^{(m+k-1)(u+\tau)+\beta/2+\eta}X^\varepsilon\bigr)
\]
by the preceding prescribed-pool test at each of the \(mk\) positions. Repeated anchors cost at most \(O_{m,k}(n^{(m+k-1)(u+\tau)})\). The power-mean lower bound for all ordered homomorphisms is
\[
\frac{|\mathcal E|^{mk}}{H_U^{mk-m}H_V^{mk-k}}.
\]
Combining the inequalities gives
\[
|\mathcal E|\ll n^{\,2u-(u-\beta/2)/(mk)+2\tau+\eta}X^\varepsilon,
\]
upon adjusting the two input tolerances. This is a fixed saving on every \(u\)-sector bounded away from zero and \(\beta=2u\), even for weights growing polynomially in \(1/L\). For probe heights \(P_h/C_h\le n^{A_h(L)}\), the same conclusion holds when \(L\max_h A_h(L)\|\lambda\|_1\|\mu\|_1\to0\), or under each exact depth comparison without an outer limit.
