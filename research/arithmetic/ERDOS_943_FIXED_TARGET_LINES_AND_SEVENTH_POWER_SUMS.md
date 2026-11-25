# Fixed-target line bounds and square-full sums at seventh powers

For a prescribed positive integer target, the line method below gives a coefficient-uniform cube-divisor bound. After gcd normalization it proves r(M^7)<<M^(2+epsilon) for every positive integer M, including noncoprime square-full summands. The unrestricted all-target improvement remains open.

## 1. Setup and conclusions

Let n>=2 and let kappa_1,kappa_2 be fixed positive integers. Count positive integer solutions
\[
 \kappa_1a^2b^3+\kappa_2c^2d^3=n,\qquad
 \gcd(\kappa_1a^2b^3,\kappa_2c^2d^3)=1,                 \tag{1}
\]
with b,d squarefree and
\[
 a\le R_1,\quad c\le R_2,\quad b\le K_1,\quad d\le K_2,
 \qquad 1\le R_i,K_i\le n.
\]
Denote their number by Z(n). All constants below are uniform in the fixed coefficient pair; they do not authorize aggregation over varying pairs. A nonempty set automatically has kappa_i<=n, gcd(kappa_1,kappa_2)=1, and every variable and coefficient in (1) is a unit modulo n. In particular
\[
 (a,c)=(b,d)=(a,d)=(c,b)=1.                            \tag{2}
\]

We first prove the two line bounds used in the unrestricted estimate and seventh-power theorem of Sections9--10.

**Lemma 1 (Fixed nonzero level on a rational line).** For any rank-two integral lattice M contained in Z^4, the number of points (a,c,b,d) in M satisfying (1), the positivity conditions, and (2), is O_epsilon(n^epsilon). The constant is independent of M and of all coefficients of a lattice basis. No bound on the coefficients of M is needed for this lemma.

**Theorem 2 (Every cube-divisor modulus).** Let q>=2 be ANY integer satisfying
\[
                              q^3\mid n.               \tag{3}
\]
No squarefreeness or lower size assumption on q is imposed.
Then
\[
 Z(n)\ll_\epsilon n^\epsilon
 \left(q+\frac{(R_1R_2K_1K_2)^{2/3}}q\right).           \tag{4}
\]
This counts representations, before choosing a graph edge representative. It applies to arbitrary fixed positive coefficients in (1); it does not require them to be squarefree.

Theorem 2 is a restricted-target result. It does not prove the full critical rectangular cycle bound, a bound for all normalized targets, or an unrestricted improvement of the prescribed-sum exponent 3/10.

The remaining fixed-target transfer issue is made precise in Sections 6--7. Heath-Brown DOES prove a lattice-frequency theorem, Lemma 6.5. Its auxiliary generator quotient is variable even when n is fixed. An explicit minimal-basis example below disproves replacing that quotient by the fixed target quotient.

## 2. A coefficient-uniform fixed-level lemma

Choose a Z-basis g,h of M, and write its four coordinates as integral linear forms
\[
 a=L_1(u,v),\quad c=L_2(u,v),\quad
 b=M_1(u,v),\quad d=M_2(u,v).
\]
Any counted point has gcd(u,v)=1: a common divisor would divide both a and c. The coordinate map is injective. If a coordinate form is identically zero there are no positive points. Define the homogeneous binary quintic
\[
 F=\kappa_1L_1^2M_1^3+\kappa_2L_2^2M_2^3.              \tag{5}
\]
If F is identically zero there are no points at the nonzero level n.

First suppose that none of the four pairs
\[
 (L_1,L_2),\quad(M_1,M_2),\quad(L_1,M_2),\quad(L_2,M_1)
                                                               \tag{6}
\]
is proportional over Q. The two degree-five forms
\[
 A=\kappa_1L_1^2M_1^3,\qquad B=\kappa_2L_2^2M_2^3
\]
have no common projective root. Each has at most two distinct roots, including a possible root at infinity. The rational map A/B:P^1_C->P^1_C therefore has degree five. The ramification contributions above zero and infinity total at least
\[
 (5-2)+(5-2)=6.
\]
If F=A+B has r distinct projective roots, its ramification contribution above -1 is 5-r. Riemann--Hurwitz gives
\[
 6+(5-r)\le 2\cdot5-2=8,
\]
so r>=3. This also covers L_1 proportional to M_1 or L_2 proportional to M_2; those coincidences increase the zero/pole ramification rather than reduce it. Only the cross-pair coincidences in (6) were excluded.

Equivalently one can use the homogeneous Jacobian of A and B, of degree eight: its factors at zeros/poles contribute at least six, and the repeated roots of F contribute 5-r, on disjoint root sets. The Jacobian is nonzero because A/B is nonconstant. This supplies the same elementary inequality without a genus computation.

Evertse's coefficient-uniform Thue--Mahler theorem applies to any integral binary form of degree five with at least three distinct complex projective roots. The statement counts all primitive integer pairs with |F(u,v)| supported on a specified set of t primes, and bounds their number by
\[
                         7^{125(2t+3)}.                \tag{7}
\]
Take the primes dividing n. The equation F(u,v)=n is a subset of that set. Thus its number of primitive solutions is at most (7) with t=omega(n), hence O_epsilon(n^epsilon). The theorem is uniform in all coefficients of F, allows reducible forms and repeated roots, and does not assume a nonzero discriminant. A form with content not dividing n has no solutions; otherwise no separate content factor is necessary. This proves the required uniformity in the nonproportional case.

It remains to treat the pairs in (6).

If L_1 and L_2 are proportional, the positive coprime pair (a,c) is fixed at most once. For that pair, b,d solve
\[
 \kappa_1a^2b^3+\kappa_2c^2d^3-nz^3=0,\qquad z=1.
\]
Heath-Brown's Theorem 1.3, with height at most n, gives n^epsilon possibilities, uniformly in the three nonzero signed coefficients. Its primitivity condition holds because z=1.

If M_1 and M_2 are proportional, the positive coprime pair (b,d) is fixed at most once. Put U=kappa_1b^3 and V=kappa_2d^3. The equation Ua^2+Vc^2=n has n^epsilon positive solutions uniformly in U,V: the algebraic integer Ua+c sqrt(-UV) has norm Un in the maximal order of Q(sqrt(-UV)); there are at most divisor-many ideals of norm Un and at most six unit generators per principal ideal of that norm. Nonemptiness gives U,V<=n, so Un<=n^2. The same argument includes the case that the field is Q(i) or Q(sqrt(-3)).

If L_1 and M_2 are proportional, the positive coprime pair (a,d) is fixed at most once. In the two-dimensional parameter space their prescribed values cut out an affine rational line, unless at most one point is left. On that affine line c and b are affine linear functions of one parameter. The target equation is
\[
 \kappa_1a_0^2b(t)^3+\kappa_2d_0^3c(t)^2=n.            \tag{8}
\]
If b is nonconstant the left side has nonzero cubic coefficient, since its other summand has degree at most two. If b is constant and c nonconstant it has nonzero quadratic coefficient. If both are constant, all four coordinates are fixed, and injectivity of the lattice parametrization permits at most one point. Thus (8) gives at most three points. Interchanging the two summands treats L_2 proportional to M_1.

These cases are exhaustive. In particular no assertion about a constant affine restriction is left implicit. Lemma 1 follows.

## 3. The fixed-q congruence and bilinear cover

For every solution (1) and any q as in (3), define the unit residue kappa modulo q by
\[
                         ab+\kappa cd\equiv0\pmod q.
\]
There are at most q choices. Put a_1=kappa_1 and a_2=kappa_2 in the following formulas, keeping kappa as the residue parameter. Heath-Brown's Lemma 4.1 derives from (1) modulo q^3 the exact constraints
\[
 a_2a\equiv\kappa^3a_1c\pmod q,\qquad
 \kappa^2a_1b+a_2d\equiv0\pmod q,                       \tag{9}
\]
\[
 2a_2ab+\kappa^3a_1cb+3\kappa a_2cd\equiv0\pmod {q^2},
                                                               \tag{10}
\]
\[
 5\kappa^2a_1a_2ab+a_2^2ad+\kappa^5a_1^2cb
                     +5\kappa^3a_1a_2cd\equiv0\pmod {q^3}.
                                                               \tag{11}
\]
Here is the direct forward derivation for arbitrary q. Equations (9) follow modulo q from kappa=-ab/(cd) and the target equation. Put
\[
 u_1=a_2a,\quad u_2=\kappa^3a_1c,\quad
 v_1=\kappa^2a_1b,\quad v_2=a_2d .
\]
These are units modulo q, with u_1=u_2 and v_2=-v_1 modulo q. Their units are also invertible modulo q^3. Choose alpha,beta in qZ/q^3Z so that
\[
 u_1/u_2=1+\alpha,\qquad -v_2/v_1=1+\beta .
\]
The target congruence is
\[
 (1+\alpha)^2=(1+\beta)^3\pmod {q^3}.
\]
It first gives 2alpha=3beta modulo q^2. Multiplication by alpha and beta gives 2alpha^2=3alpha*beta and 3beta^2=2alpha*beta modulo q^3. Doubling the original expanded identity therefore yields
\[
                  4\alpha-6\beta-\alpha\beta=0\pmod {q^3}.
\]
Multiplication by u_2v_1 gives
\[
 5u_1v_1+u_1v_2+u_2v_1+5u_2v_2=0\pmod {q^3},
\]
which is (11). Likewise 2alpha=3beta modulo q^2 gives (10). There is no division by2 or3 and no squarefree-modulus step. The converse part of Heath-Brown's Lemma4.1 is not used.

This derivation requires neither a third varying square-full summand nor the critical-range assumptions from his Section3.

Let V=R_1R_2K_1K_2. In the proof of his Lemma 4.2, before the use of his equation (3.7), the dual congruence lattice Lambda in R^4 has determinant 2q^-6. The product-coordinate scaling
\[
 E=\operatorname{Diag}(R_1K_1,R_1K_2,R_2K_1,R_2K_2)
\]
gives det(E Lambda)=2V^2q^-6. Lemma 2.8, with k=4, therefore supplies
\[
                       O(1+V^{2/3}/q^2)               \tag{12}
\]
bilinear relations
\[
          t_1ab+t_2ad+t_3cb+t_4cd=0                    \tag{13}
\]
covering all actual points with this fixed (q,kappa). The plus one in (12) is retained. Dropping it is not valid in the present shapes.

The relation vectors may be chosen integral and polynomially bounded in n: multiplying an element of Lambda by q^3 clears denominators, and Lemma 2.8 bounds the scaled vectors by an absolute constant times 1+V^{2/3}/q^2. Since all caps are between one and n, q<=n^(1/3), and the entries of E^-1 are at most one, one may use the crude bound
\[
             \max|t_i|\ll q^3+V^{2/3}q\ll n^3.        \tag{14}
\]
Only a fixed polynomial bound is needed. Inverses of a_2 modulo q^3 exist by (1). There is no dependence on the numerical sizes of the positive coefficients.

If Delta=t_1t_4-t_2t_3 is nonzero, the primitive kernel pair (b,d) gives
\[
 hb=-t_2a-t_4c,\qquad hd=t_1a+t_3c,\qquad h\mid\Delta,
                                                               \tag{15}
\]
with a possible signed nonzero h; allowing either sign changes a constant. The divisibility follows from gcd(a,c)=1. Thus at most n^epsilon rational projective lines arise per relation (13), since Delta is polynomially bounded.

For each rational projective line (15), take ALL its integral points. They form a rank-two integral lattice in Z^4. Lemma1 bounds its actual fixed-level points by n^epsilon. No further congruence refinement is necessary: retaining extra points on the rational line only enlarges the upper bound. In particular neither Lemma2.6 nor any squarefree-modulus factorization is imported here.

Consequently the nondegenerate relations contribute
\[
                   O_\epsilon(n^\epsilon(1+V^{2/3}/q^2))
                                                               \tag{16}
\]
actual level-n points for each kappa. The bound has no dependence on the eventual basis coefficient heights, by Lemma 1.

## 4. Degenerate bilinear relations and proof of Theorem 2

If Delta=0, a nonzero relation (13) factors into a rational linear form in (a,c) times one in (b,d). It fixes one of those positive primitive pairs at most once.

If (a,c) is fixed, apply the signed primitive cubic theorem directly to
\[
 \kappa_1a^2b^3+\kappa_2c^2d^3-nz^3=0,\qquad z=1 .
\]
There are n^epsilon possibilities, uniformly in all nonzero coefficients. If (b,d) is fixed, the imaginary-quadratic norm argument from Lemma1 gives n^epsilon positive root pairs. Neither argument requires a modulus size assumption.

The degenerate part therefore satisfies (16) for every q in (3). Summing over at most q residues kappa proves (4). No quadratic-congruence root count, primitive-fraction packing bound or squarefreeness of q is needed.

All assertions concern actual primitive points. The squarefreeness of b,d remains part of the original point set, while the cross-coprimalities used above follow from (1).

## 5. Exact critical specialization

Use the current critical shaped caps
\[
 R_1=\frac{2n^{11/51}}{(\kappa_1^2\kappa_2)^{1/3}},\quad
 R_2=\frac{2n^{11/51}}{(\kappa_1\kappa_2^2)^{1/3}},
\]
\[
 K_1=\frac{2n^{32/153}}{(\kappa_1^4\kappa_2^2)^{1/9}},\quad
 K_2=\frac{2n^{32/153}}{(\kappa_1^2\kappa_2^4)^{1/9}}.
\]
A cutoff below one makes the set empty. Exactly
\[
 V=\frac{16n^{130/153}}{(\kappa_1\kappa_2)^{5/3}},\qquad
 V^{2/3}=\frac{16^{2/3}n^{260/459}}
                            {(\kappa_1\kappa_2)^{10/9}}. \tag{17}
\]
Therefore the divisibility (3) and
\[
 \frac{n^{125/459}}{(\kappa_1\kappa_2)^{10/9}}
                        \le q\le n^{135/459}            \tag{18}
\]
imply Z(n)<<n^{5/17+epsilon}. Fixed absolute factors in these endpoints only alter the implicit constant. A coefficient-independent sufficient interval is
\[
                         n^{125/459}\le q\le n^{135/459}.
                                                               \tag{19}
\]
At the balanced optimum q of size V^{1/3}, the bound is n^epsilon V^{1/3}, which has coefficient-one exponent 130/459<5/17.

For example, every prime power n=p^e with e>=46 admits an integer k with
\[
                    125e/459\le k\le135e/459.
\]
Indeed the interval has length10e/459>=1. Take q=p^k. Then 3k<=e, because135/459<1/3, and (19) holds. Thus Theorem2 bounds the CRITICAL rectangular representation count at every such normalized prime power by n^{5/17+epsilon}, uniformly in p and e. This is a representation count, not a six-cycle estimate. No squarefree-q hypothesis may be reintroduced in this example.

This is a bound for the prescribed NORMALIZED target in the indicated cube-divisor sector. It does not assert that every divisor normalization of an arbitrary original square-full sum lies in that sector. Nor does a representation bound on this sector imply the full six-cycle criterion there. No new cycle exponent is being inferred from (4).

## 6. What Section 6 really counts: a concrete obstruction to fixing the generator quotient

Heath-Brown's Lemma 6.1 is a per-lattice count of square values of a binary quintic. His Lemmas 6.2--6.5 are a genuine frequency argument, not merely completion estimates. The proof of Lemma 6.2, printed page 28, writes for an AUXILIARY lattice vector g=(x_1,x_2,y_1,y_2)
\[
 18(a_1x_1^2y_1^3+a_2x_2^2y_2^3)=f q^3 r^2,\quad
 r=\gcd(y_1,y_2,q),                                    \tag{20}
\]
and sums over the nonzero integer f within a height bound. It uses Theorem 1.3 only AFTER counting x_1,x_2,r,f. Lemma 6.5 then combines counts of two independent basis vectors, together with the line-cover cap, by Cauchy--Schwarz.

Fixing an actual target n makes the left side 18n for an actual point. It does NOT do so for a generator of the containing lattice. This remains false for a successive-minimum basis, not just for an arbitrary large basis.

Here is an exact example. Take q=5, kappa=2, a_1=a_2=1 and
\[
 v=(11,2,1,1),\quad h=(125,125,125,0),\quad
 M=\mathbb Zv+\mathbb Zh,\quad
 D=\operatorname{Diag}(11,2,1,1).
\]
The actual point v satisfies
\[
 F(v)=11^2+2^2=125,\qquad n=q^3,
\]
and all primitive and squarefree conditions (the two kernels are both one; a graph cycle is not being claimed). The congruences (9)--(11) hold at v: the nontrivial expressions in (10),(11) are 50 and375. They hold at every point of M, since h is divisible by q^3 and the congruences are homogeneous linear or quadratic. Both coordinate-row maps have full rank: their determinants on (v,h) are1125 and-125. Thus every nonzero lattice vector has a nonzero root pair and a nonzero kernel pair, as required of the retained lattices in Lemma4.2.

The first normalized successive minimum is one, attained by v. Any vector uh+kv with u!=0 has last coordinates (125u+k,k), so its D-normalized sup norm is at least125|u|/2. For |u|>=2 it is at least125. For u=1 the minimum integer value is63, attained by k=-62 or-63. The first two coordinates do not increase those two norms. Thus
\[
 \lambda_1=1,\qquad\lambda_2=63,
\]
and v together with either of the following vectors is a successive-minimum basis:
\[
 h-62v=(-557,1,63,-62),\qquad
 h-63v=(-568,-1,62,-63).
\]
Their quintic values are respectively
\[
 F(h-62v)=77\,576\,593\,375=125\cdot620\,612\,747,
\]
\[
 F(h-63v)=76\,890\,082\,625=125\cdot615\,120\,661.          \tag{21}
\]
Their r in (20) is one. Neither generator quotient is the actual target quotient. These statements are exact integer identities, not a finite-survey asymptotic.

This example does not claim that a particular implementation of the finite covering algorithm selects precisely M. It disproves the purported inference from all the retained lattice/congruence properties, an actual level-n point, and successive-minimum normalization to a fixed generator value. A theorem about the algorithm's special selection would need an additional proved property.

## 7. One explicit missing transfer lemma

The fixed-level gain in Lemma1 removes the need for a square-sieve factor on each rational line. What is still missing for a general fixed target is an economical bound for the number of OCCUPIED lines produced by the congruence construction, retaining the target condition.

More precisely, fix q, the coefficients, the four shaped caps, and one finite line/lattice cover built as in Section3. For each (q,kappa) let L_n(q,kappa) be the subcollection containing at least one positive primitive point of (1). Assign a point to just one covering line if needed; do not count a point once per possible cover. The missing fixed-target frequency input is, for example,
\[
 \sum_{\kappa\in(\mathbb Z/q\mathbb Z)^\times}
       \#L_n(q,\kappa)+
       \#\{\hbox{actual points on degenerate relations}\}
                  \ll_\epsilon n^{5/17+\epsilon}.       \tag{FT}
\]
It is enough to establish this with the explicit congruence cover, not for all rational lines. The coefficient and shape hypotheses are (1), q>=2 with q^3 dividing n, and the critical cutoffs of Section5 (or another explicitly stated set of caps). Lemma1 would then prove the representation target in that rectangle. Theorem2 proves (FT) in the interval (18), by counting all lines.

Outside that interval the existing unconditional cover only supplies
\[
       \ll n^\epsilon(q+V^{2/3}/q).
\]
For q=n^{1/3}, its unavoidable displayed baseline q is n^{1/3}, larger than n^{5/17} by n^{2/51}. For targets with no nontrivial cube divisor, the only such q is one; the congruence construction has no q^3 leverage and the scale V^{2/3} is much too large. A treatment of arbitrary targets therefore also needs a different source of modulus if this route is to be used.

The statement (FT) is an unresolved frequency input, not a claimed bound. As an unweighted occupied-cover count it is close to a reformulation of the original point frequency, because each occupied line contributes only n^epsilon actual points. Its potential arithmetic content must come from new control on its basis/quotient distribution, not from merely renaming points as lines.

In particular, the natural attempt to apply Lemma6.5 while replacing f in (20) by18n/q^3 is invalid by (21). A valid refinement would have to count the joint incidence
\[
 \{(g,h):\text{congruence lattice conditions},\
          \exists(u,v)\ {\rm primitive},\
          F_{g,h}(u,v)=n\},
\]
including BOTH generator quotients, their possible zero values, and their actual size ranges. The source's global generator bounds do not condition on that exact-level incidence. This identifies a concrete additional lemma rather than assuming independence between a fixed target and lattice shape.

The path-root gcd source already counts the coefficient triples (r_1,r_2,ell) before invoking the cubic theorem. No improvement of those frequencies follows just from the per-line or primitive-cubic estimates above. The coupled weights alpha^4 gamma and alpha^3/sqrt(e), and the existing lcm/product-ratio sectors, are neither reproved nor enlarged here.

## 8. Classical inputs and source scope

Heath-Brown, *Counting Square-full Solutions to x+y=z*, arXiv2601.07817v1:
https://arxiv.org/abs/2601.07817v1 .
The exact mathematical inputs are Theorem1.3, the coefficient-uniform primitive signed diagonal-cubic bound, and Lemma2.8, the projective lattice-slicing theorem. The dual determinant and rational-line extraction follow the proof of Lemma4.2, printed pp20--22, retaining its first plus-one bound before the critical-range simplification. The forward lifting formulas of Lemma4.1 are derived directly in Section3 for every integer q. Its converse, the squarefree-q quadratic refinements, and the square-value estimate of Lemma6.1 are not needed.

The discussion of the unresolved frequency step refers precisely to Lemmas6.2--6.5, printed pp27--30: generator values vary, and two-generator frequency is a separate argument from counting values on one line. The global triple-count theorem is not imported as a fixed-target estimate.

Evertse, *Upper bounds for the numbers of solutions of certain Diophantine equations*, Seminaire de Theorie des Nombres, Paris1983--84, Birkhauser1985, pp71--90:
https://ir.cwi.nl/pub/2463/2463D.pdf .
The exact input is Theorem6, printed p84, with its supplied derivation on p85. It applies to integral binary forms of degree at least3 with at least three distinct complex projective roots, including reducible forms and repeated roots. The arguments are primitive integer pairs; for a specified t-element prime set it gives at most 7^(degree^3*(2t+3)) solutions with values supported on that set. In degree5 this is (7), independent of the coefficients. Its supplied proof imports the S-unit statement called Theorem5 in that source.

The existing coefficient-uniform Evertse1997 result for irreducible forms is a different source. No irreducibility assumption is silently substituted for the three-distinct-root theorem used here.

The critical-shaped and coupled-path sources motivate the comparison in Sections5--7. They are not premises of Theorems3--4: all dyadic resource and gcd-normalization arguments for those two theorems are written below. No checker or finite survey is a mathematical input.

## 9. Removing the rectangle: a coefficient-uniform primitive bound

Write R_kappa(n) for the number of ALL positive solutions (1), with b,d squarefree and coprime summands, without prescribed dyadic caps.

**Theorem 3 (Unrestricted fixed-target cube-divisor bound).** For every n>=2, every fixed pair of positive integers kappa_1,kappa_2, every integer q>=1 with q^3 dividing n, and every epsilon>0,
\[
 R_\kappa(n)\ll_\epsilon n^\epsilon
 \left(q+\frac{n^{1/2}}
 {(\kappa_1\kappa_2)^{1/4}q^{3/4}}\right).              \tag{22}
\]
The constants are uniform in n,q and the coefficient pair. This is not an aggregate over coefficient pairs.

First take q>=2. Partition the four variables into dyadic intervals, clipping their upper caps at n if necessary. Every nonempty box has upper caps R_1,R_2,K_1,K_2 between one and n, and its variables exceed half of their respective caps. Hence
\[
 \kappa_1 R_1^2K_1^3\le32n,\qquad
 \kappa_2 R_2^2K_2^3\le32n.                            \tag{23}
\]
Put a=R_1R_2, b=K_1K_2 and c=\kappa_1\kappa_2. Multiplication gives
\[
                       a^2b^3\le1024n^2/c.           \tag{24}
\]
Theorem2 bounds the box by n^epsilon(q+(ab)^{2/3}/q). Independently, fixing the two kernels gives n^epsilon root pairs by the same imaginary-quadratic norm argument used in Lemma1. There are at most b kernel pairs. Thus the box has at most
\[
 n^\epsilon\min\{b,\ q+C n^{2/3}c^{-1/3}q^{-1}b^{-1/3}\},
                                                               \tag{25}
\]
where C is an absolute constant. This uses (24):
\[
 (ab)^{2/3}\le1024^{1/3}n^{2/3}c^{-1/3}b^{-1/3}.
\]
For any b,L>0,
\[
                  \min(b,Lb^{-1/3})\le L^{3/4}.
\]
Indeed compare b with L^{3/4}. Also min(b,q+x)<=q+min(b,x). Applying these two elementary inequalities to (25) proves the bound (22) for a box. There are O(log(n)^4) boxes, absorbed by beginning with a smaller epsilon. No critical-locus or six-cycle estimate has been used.

For q=1 an elementary square-full count suffices, and is slightly stronger. Each first summand kappa_1*a^2*b^3 determines its second summand, and a positive square-full integer has exactly one representation c^2*d^3 with d squarefree. Therefore
\[
 R_\kappa(n)\le
 \sum_{b\le(n/\kappa_1)^{1/3}}
       \left\lfloor\sqrt{n/(\kappa_1b^3)}\right\rfloor
 \le\zeta(3/2)\sqrt{n/\kappa_1}.
\]
Swap the summands and use
\(\kappa_{\max}\ge\sqrt{\kappa_1\kappa_2}\) to obtain
\[
 R_\kappa(n)\ll
       n^{1/2}/(\kappa_1\kappa_2)^{1/4}.
\]
This proves (22) also for q=1. It uses no old three-tenths theorem.

For n=p^{7j}, taking q=p^{2j} in (22) gives the primitive bound n^{2/7+epsilon}. The next theorem includes all noncoprime summands and arbitrary seventh-power bases.

## 10. All seventh-power prescribed targets

Let r(N) count ordered positive square-full integers m_1,m_2 with m_1+m_2=N. No coprimality condition is imposed on them.

**Theorem 4 (All seventh powers).** Uniformly for every positive integer M,
\[
                  r(M^7)\ll_\epsilon M^{2+\epsilon}.
                                                               \tag{26}
\]
Equivalently, \(r(N)\ll_\epsilon N^{2/7+\epsilon}\) on the entire prescribed-target family N=M^7. The base M need not be squarefree, a prime, or a prime power. This is not a bound for arbitrary N.

Here is the full noncoprime normalization. Given an original pair, put h=gcd(m_1,m_2) and n=N/h. The integer h is square-full, since each common prime has valuation at least two in each summand. The quotients v_i=m_i/h are coprime. Let kappa_i be the product of primes having exponent exactly one in v_i. Removing those primes leaves a square-full positive integer. Consequently
\[
 v_1=\kappa_1a^2b^3,\qquad
 v_2=\kappa_2c^2d^3,
\]
with unique squarefree b,d and positive a,c, and this is a representation counted by R_kappa(n).

Each exponent-one prime in v_i divides h: its original valuation in m_i was at least two. Thus kappa_i divides rad(h); kappa_1,kappa_2 are squarefree and coprime, since v_1,v_2 are coprime. This normalization is injective for its specified h,kappa_1,kappa_2. We only need its direction from actual original pairs to normalized tuples; some arbitrary normalized tuples may fail to lift to square-full original summands. There are at most tau(N)^3 normalization choices.

It remains to choose an effective cube-divisor modulus for EVERY h dividing M^7. Write
\[
 M=\prod_p p^{m_p},\qquad h=\prod_p p^{t_p},\quad
             0\le t_p\le7m_p,
\]
and set
\[
 j_p=\min\left(2m_p,
          \left\lfloor(7m_p-t_p)/3\right\rfloor\right),
       \qquad q=\prod_p p^{j_p}.                       \tag{27}
\]
Then q is an integer, q^3 divides n=M^7/h, and q<=M^2. Most importantly,
\[
                         q\ge M^2/h^{1/2}.            \tag{28}
\]
To check (28) prime by prime, it suffices to show j>=2m-t/2. If m=1,t=0 or1, the chosen j is2, so the inequality is direct. In every other case 2m+t>=4, and
\[
 \left\lfloor(7m-t)/3\right\rfloor
       \ge(7m-t-2)/3\ge2m-t/2.
\]
The other entry 2m in the minimum is also at least2m-t/2. This proves the desired valuation inequality in all cases, including cases where its right side is negative.

Apply Theorem3 with this q, including its elementary q=1 branch. Ignoring the favorable coefficient factor is legitimate because kappa_i>=1. Equations (27)--(28) give
\[
\begin{split}
 q+\frac{n^{1/2}}{q^{3/4}}
 &\le M^2+
 \frac{M^{7/2}h^{-1/2}}
      {(M^2h^{-1/2})^{3/4}}\\
 &=M^2+M^2h^{-1/8}
 \le2M^2.
\end{split}                                                   \tag{29}
\]
Hence every nonempty normalized coefficient pair contributes at most
\(C_\epsilon N^\epsilon M^2\). The case n=1 has no pair of positive quotients and contributes zero. Summing over at most tau(N)^3 choices and using divisor bounds proves (26), after renaming epsilon. M=1 is empty and harmless.

The proof uses only the fixed-target line theorem, its coefficient-uniform Thue--Mahler/cubic/norm inputs, the elementary dyadic resource inequality and the displayed gcd normalization. It does not use the previously proved global three-tenths estimate, the critical rectangular criterion, any cycle hypothesis, or a conjectural auxiliary-cover frequency estimate.

The unsolved all-target frequency problem is untouched. In particular no assertion that a general positive integer has a suitably sized cube divisor is made. The exponent improvement in (26) is restricted to actual seventh-power targets, with noncoprime positive square-full summands included.
