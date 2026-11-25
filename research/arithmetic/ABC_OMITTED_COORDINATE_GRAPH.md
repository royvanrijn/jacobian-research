# Omitted-coordinate graphs in the \(abc\) decomposition

The powerful-number graph applies to the dyadic boxes in the \(abc\)
exceptional-set problem. The transfer uses the modulus given by one complete
summand and a common power among selected coordinates of the other two
summands. The box comparison below gives a strict gain over the four methods
used in the first version of Browning, Lichtman and Teräväinen and the box
estimates used in the later joint paper with Bernert.
[Powerful-number graph, codegree estimate](POWERFUL_SUM_GRAPHS.md#codegrees-for-any-omitted-exponent) · [Browning–Lichtman–Teräväinen, §3](https://arxiv.org/html/2410.12234v1#S3) · [Bernert–Browning–Lichtman–Teräväinen, §§2–4](https://arxiv.org/html/2410.12234v2#S2)

## Box bound

Fix \(d\) and three positive integer coefficients. On side \(t\), write the
complete value as \(q_t\prod_{j=1}^d u_{t,j}^{\,j}\), where the coordinate
\(u_{t,j}\) lies in a dyadic interval of size \(U_{t,j}\). Count tuples
whose complete values form a primitive positive triple \(a+b=c\le X\).
Choose a side \(\ell\) to fix, and write \(s,t\) for the two remaining
sides. Define
\[
F_\ell=\prod_j U_{\ell,j},
\qquad
N_\ell=q_\ell\prod_j U_{\ell,j}^{\,j}.
\]
Let \(g\ge 2\) divide every exponent in chosen nonempty sets \(I,J\)
of coordinates on sides \(s,t\). Put
\[
M_s=\prod_{j\notin I}U_{s,j},\qquad
M_t=\prod_{j\notin J}U_{t,j},\qquad
H=\prod_{j\in I}U_{s,j}^{\,j/g}\prod_{j\in J}U_{t,j}^{\,j/g}.
\]
For each \(\varepsilon>0\) the number of tuples is
\[
 \ll_{d,\varepsilon}X^\varepsilon F_\ell
 \min\!\left\{M_sM_t,\,
 M_s+M_t+(1+H/\sqrt{N_\ell})
              \min(M_s\sqrt{M_t},M_t\sqrt{M_s})\right\}.
\]
Absolute constants in dyadic interval errors depend only on \(d\). The
estimate requires pairwise coprimality between the complete values, which
follows from primitivity and the additive equation. The anatomic
decomposition can include an exceptional coefficient and coordinates that
share prime factors within a complete value.

### Proof

Fix the tuple on side \(\ell\) and call its complete value \(n\). On each
free side, choose a single tuple for each complete value represented in its
box. Replacing the free tuples by these choices preserves the complete
values and their primitivity. Each complete value has
\(O_{d,\varepsilon}(X^\varepsilon)\) tuples by the divisor bound, so the
number of original tuples is at most \(O_{d,\varepsilon}(X^\varepsilon)\)
times the resulting count, with a choice of smaller intermediate losses.
Swapping signs if needed, the two free sides give
\[
 P A^g+Q B^g=n.
\]
The signed products \(P,Q\) retain the coordinates outside \(I,J\) and
include the respective fixed coefficients. The positive integers \(A,B\)
are the \(g\)-th roots of the omitted products. The signs of \(P\) and
\(Q\) are fixed throughout the slice. All four factors are units modulo
\(n\) and \(\gcd(A,B)=1\).

Make a simple bipartite graph whose vertices are the values of \(P\) and
\(Q\). Adjacency requires a tuple in the slice. There are
\(O_d(M_s)\), \(O_d(M_t)\) vertices on the two sides. For fixed
\(P,Q\), the uniform binary-form bound proved below gives
\(O_{d,\varepsilon}(X^\varepsilon)\) positive pairs \((A,B)\).
Each such pair fixes the complete values and hence our chosen tuples.
The same binary-form estimate for coprime \(A,B\) also follows from
[Browning–Lichtman–Teräväinen, Lemma 3.4](https://arxiv.org/html/2410.12234v1#S3.SS3).

To bound a codegree, fix two distinct vertices \(P_1,P_2\) and, for any
common neighbor \(Q\), choose witnesses \(A_1,A_2,B_1,B_2\). Write
\(A_i=e\alpha_i\) and \(B_i=f\beta_i\), dividing each pair by its
greatest common divisor. The fraction
\[
 u/v=(\alpha_1\beta_2)/(\alpha_2\beta_1)
\]
is reduced. Indeed \(\alpha_1,\alpha_2\) are coprime, as are
\(\beta_1,\beta_2\) and each \(\alpha_i,\beta_i\). Eliminating \(Q\)
gives
\[
 (u/v)^g\equiv P_2/P_1\pmod n,\qquad
 e^g(P_1u^g-P_2v^g)=n(\beta_2^g-\beta_1^g).
\]
The number of roots of the congruence is \(O_g((2g)^{\omega(n)})\)
by the prime-power unit groups and the Chinese remainder theorem. This
is \(O_{g,\varepsilon}(X^\varepsilon)\), including \(n=1\) with its
single residue class.

Both entries of a reduced fraction are \(O_d(H)\). Distinct fractions
in a fixed residue class have \(uv'-u'v\) a nonzero multiple of \(n\).
Their vectors in the first quadrant therefore have angular separation
\(\gg_d n/H^2\). There are \(O_d(1+H^2/n)\) fractions per residue
class. For each fraction the four reduced factors have
\(O_{d,\varepsilon}(X^\varepsilon)\) possibilities by the divisor
bound. The coefficient of \(e^g\) is nonzero. Otherwise the identity
forces \(\beta_1=\beta_2=1\), hence \(P_1A_1^g=P_2A_2^g\). Our fixed
choice of tuple for each complete value would then give \(P_1=P_2\).
We recover \(e\) from the displayed identity and then \(Qf^g\) from one
of the two equations. Since \(Q\) divides that nonzero recovered integer,
the divisor bound gives \(O_{d,\varepsilon}(X^\varepsilon)\) choices of \(Q\).
The recovered integer has absolute value at most \(X\) because it
divides each complete value on the shared side.

Interchanging the two sides gives the same codegree bound
\(K\ll_{d,\varepsilon}X^\varepsilon(1+H^2/n)\) for the other vertex
class. For class sizes \(m_s,m_t\), Cauchy–Schwarz applied to degrees on
both orientations gives at most
\[
O\big(m_s+m_t+\sqrt K\min(m_s\sqrt{m_t},m_t\sqrt{m_s})\big)
\]
edges. The direct count gives at most \(m_sm_t\) edges. Account for
multiplicity, sum over the \(O_d(F_\ell)\) fixed tuples, and use
\(n\gg_d N_\ell\). Choosing the intermediate \(\varepsilon\) losses
smaller than the final one yields the claimed bound.

### Uniform binary-form bound

For fixed \(g\ge2\) and \(\varepsilon>0\), suppose \(P,Q\) are nonzero
integers and \(A,B,n\) are positive integers, all with absolute values
at most \(X\). For given \(P,Q,n\), the number of pairs with
\(PA^g+QB^g=n\) is \(O_{g,\varepsilon}(X^\varepsilon)\).
To see this, use the distinct roots of
\(T^g-(-1)^{g+1}QP^{g-1}\). The product of \(PA+\theta B\) over all
these roots \(\theta\) is \(P^{g-1}n\). If all the roots are rational,
they are integers, and any two distinct roots give two signed divisors of
\(P^{g-1}n\) that determine \(A,B\).

Otherwise choose a root \(\theta\) of degree \(h\ge2\). The integer norm
of \(\alpha=PA+\theta B\) in \(\mathbb Q(\theta)\) divides \(P^{g-1}n\)
because the complementary roots form a monic integer factor. The prime
ideal factorization of \((\alpha)\) shows that there are at most
\(\tau(|P^{g-1}n|)^h\) possibilities for this ideal. Every conjugate
of \(\alpha\) has absolute value \(O_g(X^2)\). Its nonzero norm is an
integer, so every conjugate also has absolute value
\(\gg_g X^{-2(h-1)}\).

For two elements generating the same ideal the ratio is a unit. The
vectors of logarithms of the absolute values of conjugates lie in a
box of side length \(O_g(1+\log X)\). Two distinct such vectors for
the same ideal have separation bounded below by a positive constant
depending only on \(g\). Indeed a unit whose log vector has sup
norm at most \(1\) has degree at most \(g\) and all conjugates bounded
by \(\exp(1)\), so there are only finitely many possibilities by
bounding the integer coefficients of its minimal polynomial. The
same finiteness argument bounds the number of units giving the zero
log vector, hence the multiplicity of each vector. Packing gives
\(O_g((1+\log X)^h)\) generators in the box. Finally, \(\alpha\)
determines \(A,B\) since \(\theta\) is irrational. The divisor bound
and the log factor give the required estimate.

## Strict improvement on a box with degree-one mass

Take all three side profiles equal to
\[
(x_1,x_2,x_3,x_4)=(2,1,20,9)/100,\qquad x_j=0\quad(j>4).
\]
Here \(x_j=\log_X U_{t,j}\), ignoring dyadic rounding, and the fixed
coefficients are \(1\). Each side has coordinate sum \(s=8/25\) and
weighted sum \(\sum_j jx_j=1\). The total coordinate mass is \(L=24/25\).
This exponent profile lies within the radical budget of the exceptional
set, for example for any fixed threshold \(\lambda>24/25\) with
\(\lambda<1\). The dyadic box may be placed at height comparable to
\(X\) through constant-factor changes.

Fix \(c\) and omit degree \(3\) on both other sides. Then
\(M_s=M_t=X^{3/25}\), \(H=X^{2/5}\), and \(N_\ell\asymp X\). The
graph exponent is
\[
s+\frac32(s-x_3)=\frac12.
\]

| Box estimate | Best exponent on this profile |
| --- | ---: |
| Fix two complete sides | \(16/25\) |
| Fourier, either version | \(27/50\) |
| Geometry of numbers, either version | \(63/100\) |
| Determinant method, first version | \(47/75\) |
| Grouped Thue bound, first version | \(14/25\) |
| Omitted-coordinate graph | \(1/2\) |

For both Fourier estimates the best grouped mass on each side is \(x_3\),
so their common exponent is \(2s-x_3/2\). The Thue estimate removes this
mass from two sides and gives \(L-2x_3\). The determinant estimate removes
\(x_j+x_k-\min(x_j/k,x_k/j)\) from \(L\) for any pair of free degrees
\(j,k\). The largest removal is \(1/3\), at \(j=k=3\). All other degree
pairs remove at most \(x_3+x_4=29/100\).

For the geometry-of-numbers estimate write \(S\) for the mass of all
omitted coordinates on the three sides and \(W\) for their weighted mass.
The saving is \(\min(S,1-W+S)\). If zero degree-three coordinates are
omitted and \(S>33/100\), all three degree-four coordinates must be
omitted, giving \(W-S\ge 81/100\). If exactly one degree-three coordinate
is omitted, then \(S\le 29/100\) without a degree-four omission, and
\(W-S\ge 67/100\) with one. If at least two degree-three coordinates
are omitted, then \(W-S\ge80/100\). The maximum saving is \(33/100\),
attained by omitting one degree-three coordinate, one degree-four
coordinate and two degree-one coordinates.
[Browning–Lichtman–Teräväinen, Propositions 3.1, 3.2, 3.5 and 3.6](https://arxiv.org/html/2410.12234v1#S3) · [Bernert–Browning–Lichtman–Teräväinen, Propositions 2.3, 3.1 and 4.1](https://arxiv.org/html/2410.12234v2#S2)

The comparison is between box bounds at one admissible exponent profile.
A global improvement requires controlling the maximum over the full
exceptional-set optimization.

## Comparison with the complete modulus-locking obstructions

Zeng gives two rational feasible points for his exponent programs.
The omitted-coordinate graph can be optimized exactly on both points.
[Zeng, equations (8.1) and (8.7)](https://hexagonmath.org/pdf/2611.00621v1)

Write \(s_\ell\) for the fixed side's total coordinate mass,
\(m_s,m_t\) for the two remaining vertex masses, and \(h\) for the
logarithmic root height \(\log_X H\). Both points have weighted
height \(1\) on each side. For the symmetric point (8.7), mass
\(5/396\) remains above degree \(6\) on each side, with residual
weighted height at most \(35/396\). Any lift to a full exponent
profile puts that entire mass at degree \(7\). Point (8.1) has no
tail. The graph bound has exponent
\[
s_\ell+\min\!\left\{
 m_s+m_t,\;
 \max\!\left(m_s,m_t,\min(m_s+m_t/2,m_t+m_s/2)+(h-1/2)_+\right)
 \right\}.
\]
It is nondecreasing in each of \(m_s,m_t\) and in the nonnegative
height penalty. A lower bound on its minimum follows by removing
the largest possible mass in degrees divisible by one integer \(g\ge2\)
on each side and setting the height penalty to zero.

At the symmetric point, the largest removable mass on each side
is \(3/22\), obtained in degrees \(2,4\) or degree \(3\). A larger
degree removes at most \(1/11\), including the option of omitting
degree \(7\). Thus \(m_s,m_t\ge13/66\). Taking degree \(3\) on both
sides gives \(h=3/11<1/2\) and attains the lower bound
\[
 \frac13+\frac32\frac{13}{66}=\frac{83}{132}.
\]
The distance from the point's \(D=7/12\) to this right-hand side
is \(1/22\). The point still satisfies every full-side
common-power graph inequality.

At the old point (8.1), the greatest removable mass is \(3/23\) on
\(a,b\) and \(53/345\) on \(c\). These three maxima may all be
realized with degrees \(2,4\). Fixing \(a\) or \(b\) gives a lower
bound for the graph exponent with fixed mass \(37/115\) and other
vertex masses \(22/115,14/69\). This choice has \(h=32/69<1/2\),
so the attained minimum is
\[
 \frac{37}{115}+\frac{22}{115}+\frac12\frac{14}{69}
   =\frac{212}{345}.
\]
Fixing \(c\) permits exponent at least
\(41/115+(3/2)(22/115)=74/115>212/345\).
The old point has \(D=403/690\) and satisfies the graph inequalities
with a minimum margin of \(7/230\).


## Collision bounds for square and cube roots drawn from sparse sets

In a fixed anatomic box, let \(A,B,C\) be the three letters, with coordinate
scales \(X^{\alpha_{T,i}}\) on letter \(T\). Write \(L_T=\sum_i\alpha_{T,i}\).
Release sets \(I_A,I_B\) of coordinates on the first two letters and freeze a set
\(O_C\) on the third. Inside \(I_A\), choose two nonempty disjoint groups
\(G_1,G_2\), with fixed positive integers \(g_1,g_2\) such that \(g_j\mid i\)
for every \(i\in G_j\). Put
\[
\begin{aligned}
 p_j&=\sum_{i\in G_j}\alpha_{A,i},
 &H_j&=\sum_{i\in G_j}(i/g_j)\alpha_{A,i},
 &r&=\sum_{i\in I_A\setminus(G_1\cup G_2)}\alpha_{A,i},\\
 s_B&=\sum_{i\in I_B}\alpha_{B,i},
 &h_B&=\sum_{i\in I_B}i\alpha_{B,i},
 &w_B&=\sum_{i\notin I_B}i\alpha_{B,i},\\
 f_A&=L_A-p_1-p_2-r,
 &f_B&=L_B-s_B,
 &f_C&=\sum_{i\in O_C}\alpha_{C,i},
 \qquad w_C=\sum_{i\in O_C}i\alpha_{C,i}.
\end{aligned}
\]
Here \(H_j\) measures the size of the grouped root, while \(p_j\) measures
the number of choices for that root. With \(x_+=\max(x,0)\), define
\[
 K=2r+p_1+p_2+\frac{(2H_1-w_B)_++(2H_2-w_B)_+}{2}.
\]
The number of primitive triples in the box has exponent at most
\[
 f_A+f_B+f_C+\max\left\{
 L_C-f_C,\frac{L_C-f_C+K+\delta}{2}
 \right\}+\varepsilon.
 \tag{1}
\]
For \(I_B\ne\varnothing\) one may take the minimum of the following values
for \(\delta\).
\[
 \begin{gathered}
 s_B,\qquad
 2\!\!\sum_{\substack{i\in I_B\\ e\nmid i}}\alpha_{B,i}\quad(e\ge2),\qquad
 (h_B-w_C)_+,\qquad
 s_B-\alpha_{B,t}+(\alpha_{B,t}-w_C)_+\quad(t\in I_B).
 \end{gathered}
\]
More generally, choose nonempty \(J\subseteq I_B\) and \(g\ge1\) dividing
every index in \(J\). A further option is
\[
 \delta_{J,g}=s_B-\sum_{i\in J}\alpha_{B,i}
           +\left(\sum_{i\in J}(i/g)\alpha_{B,i}-w_C\right)_+ .
\]
The whole-value and single-coordinate terms are the cases \((J,g)=(I_B,1)\)
and \((J,g)=(\{t\},t)\), respectively.

All exponent groups are fixed, and the constant in the bound may depend
on them, the truncation \(d\), and \(\varepsilon>0\). If only the first six
coordinates are retained individually, a lower bound for \(w_B\) from the
remaining indices may be used. This increases the right-hand side. The
formula extends the collision bound of Zeng by allowing each of the two
factors in its twisted energy to be a grouped root.
[Zeng, Lemma 3.3 and Proposition 4.1](https://hexagonmath.org/pdf/2611.00621v1)

### Multiplicative energy of sparse sets

Let \(\mathcal U,\mathcal V\) be sets of positive integers at most \(U,V\),
respectively, all units modulo \(M\). Assume \(U,V,M\) are bounded by fixed
powers of \(X\). For a unit \(\gamma\bmod M\), the number of solutions of
\[
 u^{g_1}v^{g_2}\equiv\gamma(u')^{g_1}(v')^{g_2}\pmod M
\]
is
\[
 \ll X^\varepsilon|\mathcal U||\mathcal V|
 \sqrt{(1+U^2/M)(1+V^2/M)}.
 \tag{2}
\]
Indeed, the squared \(\ell^2\) norm of the ratio counts
\(R_{\mathcal U}(k)=\#\{(u,u')\in\mathcal U^2:u/u'\equiv k\bmod M\}\)
counts quadruples with equal ratios. Fix the two factors on one side of the
cross-multiplied congruence, at a cost of \(|\mathcal U|^2\) choices. The product
on the opposite side has \(O(1+U^2/M)\) possible integer values, and each has
\(O_\varepsilon(X^\varepsilon)\) factorizations by the divisor bound.
The same argument controls \(R_{\mathcal V}\). Pushing the ratio counts
forward by \(k\mapsto k^{g_j}\) increases the squared norm by at most the
kernel size. This is \(O_{g_j,\varepsilon}(X^\varepsilon)\), since the
prime-power unit groups have at most two cyclic factors and
\(g_j^{2\omega(M)}=O_{g_j,\varepsilon}(X^\varepsilon)\). Cauchy--Schwarz
between the two pushed-forward ratio counts proves (2), uniformly in
\(\gamma\). The argument includes \(M=1\).

### From energy to the box exponent

Fix the frozen coordinates and write the additive equation, with fixed signs, as
\(\pm M_Au\pm M_Bv=\pm M_C w\) for positive monomial values
\(u,v,w\). Each value represents only divisor-many coordinate tuples, so it
suffices to count distinct value triples up to \(X^\varepsilon\). In a
nonempty fibre, the coefficients \(M_A,M_B,M_C\) are pairwise coprime.
Primitivity makes every coordinate of \(u\) a unit modulo \(M_B\)
and every coordinate of \(v\) a unit modulo \(M_C\).

For each value of \(w\), let \(n_w\) count the permitted pairs \((u,v)\), and
let \(E\) count ordered pairs of distinct permitted pairs with the same \(w\).
Writing \(N=\sum_w n_w\) and \(W=\#\{w:n_w>0\}\), Cauchy--Schwarz gives
\(N^2\le W(N+E)\), hence \(N\ll W+\sqrt{WE}\). Here \(W\ll X^{L_C-f_C}\).

An off-diagonal collision gives \(u\equiv u'\pmod {M_B}\) with \(u\ne u'\).
Freeze both copies of the coordinates in \(I_A\setminus(G_1\cup G_2)\).
This costs \(O_d(X^{2r})\) and leaves the congruence in (2). The \(j\)-th root
set consists of products \(\prod_{i\in G_j}x_i^{i/g_j}\). It has
\(O_d(X^{p_j})\) values at most \(O_d(X^{H_j})\). Every choice of grouped
roots and frozen coordinates determines \(u,u'\), so (2) bounds their number
by \(O_{d,\varepsilon}(X^{K+\varepsilon})\).

For a fixed \(u,u'\), the integer difference \(v-v'\) is fixed and nonzero.
One can count all possible \(v\), or freeze both copies of the coordinates
of \(v\) not divisible in degree by \(e\) and apply the uniform binary-form
bound above to the remaining \(e\)-th powers. These are the first two
choices for \(\delta\). One can also reduce the original additive equation
modulo \(M_C\). It determines the residue of \(v\) because \(M_B\) is
invertible there, giving the whole-value bound from \(v\ll_d X^{h_B}\).
For the grouped residue bound, freeze the coordinates of \(v\) outside
\(J\). The prescribed residue of \(v\bmod M_C\) becomes a unit \(g\)-th
power congruence in the remaining grouped root. There are at most
\(g^{2\omega(M_C)}\) admissible residue classes. The root lies in an
interval of length \(O_d(X^{\sum_{i\in J}(i/g)\alpha_{B,i}})\), and its
value determines the corresponding factor of \(v\). This gives
\(\delta_{J,g}\) after restoring frozen choices. All these estimates
count the same conditional set of \(v\)
values, so the minimum is valid. The pair \(u,v\) determines \(w\), hence
\(E\ll_{d,\varepsilon}X^{K+\delta+\varepsilon}\).
Restoring the frozen choices and the divisor losses proves (1).
[Zeng, Corollary 3.2 and Proposition 4.1](https://hexagonmath.org/pdf/2611.00621v1)

### The symmetric \(7/12\) profile

Use the symmetric profile in Zeng's equation (8.7). In units of \(1/396\),
the first six coordinate masses of each letter are \((19,18,54,36,0,0)\),
and its total and tail masses are \(132\) and \(5\). Use
\(I_A=\{1,\ldots,6\}\), \(G_1=\{2,4\}\), \(G_2=\{3\}\), \(g_1=2\),
\(g_2=3\), \(I_B=\{3\}\), and \(O_C=\varnothing\).
Both root sets have cardinality exponent \(54/396\). The square root
\(x_2x_4^2\) has size exponent \(90/396\), and the cube root \(x_3\)
has size exponent \(54/396\). The frozen \(B\)-modulus has exponent at
least \(234/396\), which exceeds both doubled root size exponents.
Consequently \(K=146/396\), \(\delta=0\), \(f_A=5/396\),
\(f_B=78/396\), and \(f_C=0\). The two terms in (1) have exponents
\(215/396\) and \(222/396\) before the arbitrary loss. Thus this box
inequality cuts the symmetric point by
\[
 \frac{7}{12}-\frac{37}{66}=\frac1{44}.
\]
[Zeng, equation (8.7)](https://hexagonmath.org/pdf/2611.00621v1)

### An asymmetric FullOC witness at smaller radical budgets

Zeng's FullOC program also has an asymmetric feasible segment on
\(243/244\le\lambda\le1\). Put \(u=244(1-\lambda)\) and divide all
entries of this table by \(1464\).

| Letter | \(\alpha_{1}\) | \(\alpha_{2}\) | \(\alpha_{3}\) | \(\alpha_{4}\) | \(\alpha_{5}\) | \(\alpha_{6}\) | Total | Fourier variable |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| \(a\) | \(80-10u\) | \(62+4u\) | \(198\) | \(142-10u\) | \(0\) | \(0\) | \(496-10u\) | \(336-12u\) |
| \(b,c\) | \(67+3u\) | \(68-2u\) | \(198\) | \(130+2u\) | \(0\) | \(0\) | \(484+2u\) | \(198\) |

The exponent coordinate is \(D=(854-2u)/1464=7/12-(1-\lambda)/3\).
The total radical budget is \(\lambda\). Realizing each tail at degree
seven gives weighted height exactly one on every letter.

The [exact witness verifier](verify_abc_fulloc_asymmetric_witness.py)
checks the root rows, all \(64^3\) geometric choices, and all \(181{,}440\)
FullOC collision instances with all 15 exponent pairs. Each affine
inequality is checked at \(u=0,1\). For each positive part or outer
maximum it uses an affine lower bound obtained by choosing the
lexicographically larger pair of constant and slope. It tests each
difference cap separately, since \(\max(d,\min_i x_i)\ge D\) follows
from \(\max(d,x_i)\ge D\) for every \(i\). This gives a lower bound on
the finite program throughout the interval. The old program cannot
give an upper exponent \(7/12-s(1-\lambda)\) with \(s>1/3\) on any left
neighborhood of \(1\). [FullOC inequalities and program, §§2, 4–6 and 8](https://hexagonmath.org/pdf/2611.00621v1)

The sparse-root bound excludes the entire segment. Choose \(A=b\),
\(B=c\), \(C=a\), \(I_A=\{1,\ldots,6\}\), \(G_1=\{2,4\}\) with
\(g_1=2\), \(G_2=\{3\}\) with \(g_2=3\), \(I_B=\{2,4\}\), and
\(O_C=\varnothing\). The difference degree \(e=2\) gives \(\delta=0\).
The two positive-part penalties vanish on the interval, so the
energy exponent is twice the mass of degree one on \(b\) plus the
masses of degrees two, three, and four. The off-diagonal term in (1)
dominates the diagonal term and gives the upper exponent
\(205/366-(1-\lambda)/6\). The excess of the proposed witness exponent
over this bound is \(17/732-(1-\lambda)/6\ge 11/488\).

### A \(4/7\) obstruction for the grouped six-index program

In units of \(1/2100\), give each letter first-six masses
\((168,0,300,200,25,0)\), total mass \(700\), and Fourier variable
\(300\). Set \(D=1200/2100=4/7\). The projected height on each
letter is \(1993+7\cdot7=2042\) in the same units. The lift with
\(5\) additional mass at degree \(15\) and \(2\) at degree \(16\)
has weighted height exactly \(2100\). The relaxed six-index
budget, pair and Fourier constraints hold by substitution.

The minima of the geometric rows and of (1) at this point are
both \(4/7\). For example, the geometric choice
\((I_a,I_b,I_c)=(\{1,3\},\{1,3\},\varnothing)\) gives maximum
\(\max(2100-936,1200)/2100=4/7\). For (1), use the degree-three
and degree-four groups on \(A\), release \(I_B=\{3\}\), and take
\(O_C=\varnothing\). The degree-three difference cap is zero,
\(K=886/2100\), and the off-diagonal term is \(4/7\).

The [exact verifier](verify_abc_grouped_fulloc_obstruction.py)
computes both minima using scaled integers. Symmetry reduces the
geometric choices to multiplicities in \(\{0,1,2,3\}\) for each
positive-mass degree. For the collision rows, it takes every pair
of disjoint subsets on \(A\), every release on \(B\), every lock on
\(C\), and every subgroup lock inside the chosen release on \(B\).
Removing a zero-mass degree from a group can only increase its
greatest common divisor and decrease its root height. The verifier
also allows empty support groups with height zero, so its search
includes a bound at least as strong as any choice involving
zero-mass degrees. For any fixed group, its greatest common
divisor gives the smallest root height. This gives a \(4/7\) lower obstruction for the finite exponent program.

The full-side common-power graph also permits this point. In the
lift, a common power can remove at most \(305/2100\) from either
free letter, attained by grouping degrees \(3\) and \(15\).
Both graph vertex masses are therefore at least \(395/2100\).
The least graph right-hand side is
\(1/3+(3/2)(395/2100)=517/840\), attained by those groups
because their combined root height is \(650/2100<1/2\).

The same lift survives projection through degree seven, and in
fact through degree fourteen. In that range the individually
retained support is unchanged. Raising the cutoff increases the
lower bound for the frozen modulus \(M_B\) and can only decrease
the collision upper bounds. The new difference degrees in this
range divide none of the positive-mass indices. The minima at
cutoff fourteen are still \(4/7\). At cutoff fifteen the collision
minimum becomes \(2327/4200\) for this lift. The cut uses
\(G_1=\{4\}\), \(G_2=\{3,15\}\), \(I_B=\{1,3,5,15\}\),
\(O_C=\{5,15\}\), and \(J=\{3\}\) for the coordinate lock.
Exposing both tail degrees gives \(2293/4200\) with
\(G_1=\{3,15\}\), \(G_2=\{4,16\}\),
\(I_B=\{1,3,5,15,16\}\), \(O_C=\{5,15,16\}\), and \(J=\{3\}\).
These latter cuts concern the indicated lift.

From the repository root run

```sh
python3 research/arithmetic/verify_abc_grouped_fulloc_obstruction.py
```


## Counting collisions by their common shift

In the notation of (1), write \(s_T=\sum_{i\in I_T}\alpha_{T,i}\) and
\(h_T=\sum_{i\in I_T}i\alpha_{T,i}\) for the free mass and weighted free
mass on \(T=A,B\). Let \(w_T\) be any lower exponent bound on its frozen
modulus. Upper bounds for \(h_T\) can be used throughout. For either free
letter, use the fixed nonzero difference cap
\[
 d_T=\min\left\{s_T,\ \min_{e\ge 2}
             2\!\!\sum_{\substack{i\in I_T\\e\nmid i}}\alpha_{T,i}\right\},
\]
where the minimum uses any chosen finite set of fixed degrees. For the
second free letter after the opposite free value has been fixed, use
the cap
\[
 \ell_T=\min\left\{s_T,\ L_C-f_C,\ (h_T-w_C)_+,\
 \min_{\varnothing\ne J\subseteq I_T}
 \left(s_T-\sum_{i\in J}\alpha_{T,i}
      +\left(\sum_{i\in J}(i/g_J)\alpha_{T,i}-w_C\right)_+\right)\right\},
 \qquad g_J=\gcd(J).
\]
A cap whose required degrees or heights are unavailable may be omitted
from the minimum. These caps count the values on one letter. Define
\[
 \sigma=\bigl(\min\{h_A-w_B,h_B-w_A,1-w_A-w_B\}\bigr)_+,
 \qquad
 \rho=\min\{d_A+d_B,\ d_A+\ell_B,\ d_B+\ell_A\}.
\]
Then the exponent of the off-diagonal collision count in a fixed fibre
is at most \(\sigma+\rho+\varepsilon\). Consequently the box exponent
is at most
\[
 f_A+f_B+f_C+
   \max\left\{L_C-f_C,\ \frac{L_C-f_C+\sigma+\rho}{2}\right\}
   +\varepsilon.
 \tag{3}
\]
For tagged moduli the same estimate holds on the tagged family. The
bounds in (1) and (3) can be minimized over the same family of boxes.

To prove the estimate, fix the frozen coordinates. Primitivity in any
nonempty fibre makes \(M_A,M_B,M_C\) pairwise coprime. Subtract the two
signed equations for a collision with common value \(w\). There is a
nonzero integer \(k\) such that
\[
 u-u'=M_B k,\qquad v-v'=\pm M_A k.
\]
The two free heights bound \(|k|\) from above, as does the complete
height bound \(M_A M_B |k|\le X\). The number of possible shifts is
therefore \(O_d(X^\sigma)\). If one of the three upper exponents inside
the minimum is strictly negative, there are no off-diagonal collisions
for sufficiently large \(X\).

Fix one shift \(k\). Each of the differences on \(A\) and \(B\) is
now fixed and nonzero. A free value determines its other copy. The
trivial cap \(s_T\) counts possible free values. Freezing both copies
of the coordinates with degrees not divisible by \(e\) and applying
the uniform binary-form bound to the two remaining \(e\)-th powers
gives the other cap in \(d_T\). Given a choice on \(A\), reduction of
the signed equation modulo \(M_C\) determines one unit residue for
\(v\), and the same holds with the letters exchanged. The whole-value
interval bound gives \((h_T-w_C)_+\). For the subgroup bound, fix the free
coordinates outside \(J\) in the chosen copy on \(T\) and solve a unit
\(g_J\)-th power congruence for the grouped root. It has \(O_{g_J,\varepsilon}
(X^\varepsilon)\) residue classes, and the weighted root height gives
the displayed interval bound. The extra term \(L_C-f_C\) counts the
possible values of \(w\), since each \(w\) fixes the second free value
by the signed equation. The caps \(d_B\) and \(\ell_B\) both apply
after choosing the value on \(A\), so they contribute
\(\min(d_B,\ell_B)\). Choosing \(B\) first gives \(d_B+\min(d_A,\ell_A)\).
This proves the paired cap \(\rho\). Once \(u,v\) are known, the signed
equation determines \(w\). Multiplying by the number of shifts and
using \(N^2\le W(N+E)\) as in (1) proves (3).
[Zeng, Corollary 3.2 and Proposition 4.1](https://hexagonmath.org/pdf/2611.00621v1)

### Releasing or freezing all higher coordinates

For a six-index model, put \(r_T=L_T-\sum_{i\le6}\alpha_{T,i}\).
One can place the entire tail in either the free or the frozen part
of a letter. When it is frozen, use
\[
 w_T\ge \sum_{i\le6,\ i\notin I_T}i\alpha_{T,i}+7r_T,
 \qquad h_T=\sum_{i\le6,\ i\in I_T}i\alpha_{T,i}.
\]
If the complete value of this letter is comparable to \(X\), the
identity \(M_Tu_T\asymp X\) also gives \(w_T\ge 1-h_T\) in the limiting
exponent. When the tail is free, include \(r_T\) in \(s_T\) and
use \(h_T\le 1-w_T\) with the \(w_T\) obtained from the early
frozen degrees. Charge all of \(r_T\) in the mass of degrees not
divisible by \(e\) when computing a grouped difference cap. A subgroup
lock can use early degrees alone, with the released tail included in
the cost of its complement. The whole-value lock uses the full
upper bound for \(h_T\). The same frozen-modulus improvement applies
when the tail of \(C\) is frozen and its complete value is comparable
to \(X\), with \(1-h_C\) obtained using just the early free degrees.
These estimates use only total tail masses. Cuts based on the
height of a given letter require that letter to be tagged in the
program. Among the original positive terms, \(a+b=c\asymp X\)
ensures both \(c\) and \(\max(a,b)\) are comparable to \(X\).
[Whole-tail modulus bound](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#using-all-higher-coordinates-of-the-largest-letter)

### The symmetric \(4/7\) profile

For the symmetric point above, release \(\{1,3\}\) on both \(A\) and
\(B\) and freeze nothing on \(C\). In units of \(1/2100\), each
free letter has \(s_T=468\), \(f_T=232\), \(h_T=1068\), and the
frozen-modulus lower bound from the six-index data is \(w_T=974\).
The cubic difference cap is \(336\) on each side. Thus
\(\sigma=94\), \(\rho\le672\), and the diagonal and off-diagonal
terms in (3) are \(1164/2100\) and \(1197/2100\). The bound is
\(399/700\), cutting the \(4/7\) point by \(1/700\) in the original
six-index variables. If either of the two free letters is comparable
to \(X\), its frozen-modulus exponent can be raised to \(1032\)
and the upper bound becomes \(292/525\).

### Direct congruence bound from two residue energies

Choose a set \(O_C\) of frozen coordinates on letter \(C\). Write
\(f_C=\sum_{i\in O_C}\alpha_{C,i}\), and let \(w_C\) be a lower bound
for the exponent of the resulting fixed coefficient \(M_C\), including
the small exceptional coefficient. Release every coordinate on letters
\(A,B\). Their full monomial products are denoted \(u,v\). For each fixed
choice on \(C\), the signed equation implies
\[
u\equiv\gamma v\pmod{M_C}
\]
for a unit \(\gamma\), because the exceptional coefficients of \(A,B\)
are coprime to \(M_C\) in any nonempty primitive box. Only unit values
of \(u\) and \(v\) contribute. Each pair of full products determines the
full value of the third letter. The divisor bound for the fixed number
of coordinates bounds the tuple multiplicity at subpower cost.

For \(T=A,B\), choose nonempty disjoint index groups
\(G_{T,1},G_{T,2}\) with fixed positive divisors \(g_{T,j}\) of their
indices. Set \(p_{T,j}=\sum_{i\in G_{T,j}}\alpha_{T,i}\) and
\(H_{T,j}=\sum_{i\in G_{T,j}}(i/g_{T,j})\alpha_{T,i}\).
Let \(E_T\) count pairs of permitted full products in the same
residue class modulo \(M_C\). Freeze both copies of the coordinates
outside the two groups. Applying (2) to the two resulting grouped root
sets, uniformly in the unit twist from the fixed coordinates, gives
\[
 E_T\ll_{d,\varepsilon}X^{K_T+\varepsilon},\qquad
 K_T=2L_T-p_{T,1}-p_{T,2}
      +\tfrac12\!\sum_{j=1}^2(2H_{T,j}-w_C)_+.
\]
Cauchy–Schwarz on the two residue multiplicity functions bounds all
congruent pairs by \(\sqrt{E_A E_B}\). Restoring the \(O_d(X^{f_C})\)
choices of the frozen coordinates gives
\[
 D\le L_A+L_B+f_C-\tfrac12\!\sum_{\substack{T=A,B\\j=1,2}}p_{T,j}
       +\tfrac14\!\sum_{\substack{T=A,B\\j=1,2}}(2H_{T,j}-w_C)_+
       +\varepsilon .
\]
The power signs do not affect the argument. All unused higher coordinates
contribute to the \(L_T\) terms, so the formula holds uniformly for their
possible distributions. If every higher coordinate on a large letter \(C\)
is frozen, one may instead use the lower bound for \(w_C\) from the known
height of that letter and the complement of \(O_C\).

At the symmetric \(4/7\) profile, freeze degree \(4\) of \(C\) and
use \(G_{T,1}=\{3\},G_{T,2}=\{4\}\) on both released letters. In units of
\(1/2100\), one has \(L_T=700,\ f_C=200,\ w_C=800\) and the root
cardinality and height exponents are \(300,200\). Both positive parts
vanish on each letter, giving \(D\le1100/2100=11/21\).
At the symmetric \(7/12\) profile use \(O_C=\{4\}\) and
\(G_{T,1}=\{3\},G_{T,2}=\{2,4\}\). In units of \(1/396\), the two
cardinality exponents are \(54,54\) and the heights are \(54,90\), with
\(L_T=132,\ f_C=36,\ w_C=144\). Thus \(K_T=174\) on each letter and
the bound is \(D\le210/396=35/66\).

### An asymmetric obstruction for the height-tagged bounds

The height-tagged projected program with whole-tail releases, grouped
collision bounds, direct residue energies and paired shift bounds has a
feasible point at \(D=2627/4796\). Its coordinate and Fourier numerators
in units of \(1/95920\) are

| Letter | \(\alpha_1\) | \(\alpha_2\) | \(\alpha_3\) | \(\alpha_4\) | \(\alpha_5\) | \(\alpha_6\) | \(L_T\) | \(F_T\) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Short \(a\) | 4508 | 5169 | 9221 | 11550 | 1266 | 0 | 31714 | 16719 |
| Tall \(b\) | 4504 | 5380 | 9552 | 11335 | 1332 | 0 | 32103 | 16715 |
| Tall \(c\) | 4504 | 5380 | 9552 | 11335 | 1332 | 0 | 32103 | 16715 |

There is no tail mass. The weighted height of \(a\) is \(95039/95920\),
and that of each tall letter is one. The six-index ordering holds
because the \(\alpha_4+2\alpha_5\) numerators are \(14082,13999,13999\).
The exact minima of the geometric, grouped collision and direct
residue-energy rows are all \(2627/4796\). The paired shift rows have
minimum \(55257/95920\). The [standard-library checker](verify_abc_tall_residue_obstruction.py)
enumerates the early and whole-tail masks, the disjoint root groups
and the possible residue locks. This is a lower obstruction for the
stated exponent program.

### Freezing a divisor of one coordinate

The direct residue argument also applies on a stratum where the degree
\(i\) coordinate \(x_{C,i}\) of the target letter has a divisor
\(m\asymp X^\theta\). Choose one such divisor per tuple, then freeze
\(m\) at cost \(X^\theta\). The complete value of \(C\) is divisible by
\(m^i\), and primitiveness makes both other letters units modulo \(m^i\).
In the residue-energy bound one can therefore take \(f_C=\theta\) and
\(w_C=i\theta\) for this lock. Each pair of full values on the released
letters determines the complete value of \(C\), so the remaining choice
of its coordinates still has only subpower multiplicity.

At the asymmetric profile above, take \(C=c\), \(i=4\) and
\(\theta=5775/95920\). Use the singleton groups \(\{3\},\{4\}\) on
each released letter. The modulus exponent \(4\theta=23100/95920\)
is at least both doubled root heights on both letters. The direct
residue bound on this stratum is then
\[
 D\le\frac{63817+5775-(20771+20887)/2}{95920}
     =\frac{48763}{95920}.
\]
The complementary stratum consists of coordinates without a divisor at
this scale. It contains prime coordinates and needs its own estimate.

At this profile the same singleton groups give \(D\le 27/50+\varepsilon\)
whenever \(x_{c,4}\) has a divisor in \([X^\ell,X^u]\), where
\(\ell=60101/1438800\) and \(u=11011/119900\). In units of \(1/95920\),
the exponent bound before the \(X^\varepsilon\) loss is
\[
42988+t+\sum_{s\in\{4610.5,4776,5667.5,5775\}}(s-t)_+,
\qquad t=95920\theta .
\]
It is convex in \(t\) and equals \(95920\cdot 27/50\) at both endpoints
\(t=95920\ell,95920u\). Dyadic subdivision of the divisor range costs
\(O(\log X)\). Since \(2\ell\le u\), when every prime factor is below
\(X^\ell\), the first partial product of prime factors to reach \(X^\ell\)
is at most \(X^{2\ell}\). A prime factor in the interval is itself a
suitable divisor. Thus for \(x_{c,4}\asymp X^{11335/95920}\), the
complement consists of \(x_{c,4}=ps\) with \(p>X^u\) prime and
\(s\ll X^{12631/479600}\), since \(11335/95920-u=12631/479600<\ell\).
For large \(X\) the cofactor is below \(X^\ell\), and every divisor
containing \(p\) is above \(X^u\).

### Higher multiplicative moments for a residue lock

Fix a modulus \(M\asymp X^w\), and let \(G=(\mathbb Z/M\mathbb Z)^\times\).
The argument also uses the one-element group if \(M=1\).
For each coordinate \(i\) take an integer set \(\mathcal X_i\) of units
modulo \(M\) in an interval contained in \([1,C_iX^{\alpha_i}]\).
The intervals and the sets need no distribution assumption.
Give this coordinate a fixed positive power \(e_i\), and put
\[
 F_i(\chi)=\sum_{x\in\mathcal X_i}\chi(x)^{e_i},
 \qquad \mathbb E_\chi f(\chi)=\frac{1}{|G|}
                   \sum_{\chi\in\widehat G}f(\chi).
\]
The scales and modulus are bounded by fixed powers of \(X\),
and the coordinate powers are fixed. Integer tuples are counted with multiplicity here.
This bounds counts of distinct monomial products because a represented
product has at least one tuple. For a fixed number of coordinates, the
divisor bound also bounds the multiplicity of each product at subpower
cost.

For an integer \(h\ge1\) fixed independently of \(X\), orthogonality
expresses \(\mathbb E_\chi |F_i(\chi)|^{2h}\) as the number of solutions to
\[
 (x_1\cdots x_h)^{e_i}\equiv (y_1\cdots y_h)^{e_i}\pmod M,\qquad
 x_j,y_j\in\mathcal X_i.
\]
There are \(O(X^{h\alpha_i})\) choices on the first side. For each,
the product on the other side lies in
\([1,O(X^{h\alpha_i})]\) and in one of
\(O_{e_i,\varepsilon}(X^\varepsilon)\) residue classes.
The latter bound follows from the size of the kernel of the
\(e_i\)-th power map on \(G\), by the prime-power unit groups and the
Chinese remainder theorem. Each possible integer product on the
second side has divisor-many factorizations since \(h\) is fixed.
Consequently
\[
 \mathbb E_\chi |F_i(\chi)|^{2h}
 \ll X^{h\alpha_i+(h\alpha_i-w)_++\varepsilon}. \tag{4}
\]
If only \(M\gg X^w\) is known, the same upper bound holds with that
lower exponent.

For a coordinate set \(J\), integer \(k\ge1\), and any unit twist
\(\gamma\), write \(E(J,k;\gamma)\) for the number of choices of
\(k\) copies of all its coordinates on each side satisfying
\[
 \prod_{\ell=1}^k\prod_{i\in J}x_{i,\ell}^{e_i}
 \equiv \gamma\prod_{\ell=1}^k\prod_{i\in J}y_{i,\ell}^{e_i}\pmod M.
\]
Orthogonality gives
\[
 E(J,k;\gamma)=
 \mathbb E_\chi\overline{\chi(\gamma)}
        \prod_{i\in J}|F_i(\chi)|^{2k}
 \le E(J,k;1).
\]
For a disjoint split \(J=J_1\sqcup J_2\), Cauchy–Schwarz therefore gives
\[
 E(J,k;\gamma)\le
 \bigl(E(J_1,2k;1)E(J_2,2k;1)\bigr)^{1/2}. \tag{5}
\]
This can be iterated to any fixed depth. For example, if each
coordinate \(i\) is placed at depth \(d_i\) in a binary split,
(4) gives the exponent
\[
 k\sum_{i\in J}\alpha_i+
 \sum_{i\in J}\bigl(k\alpha_i-w/2^{d_i}\bigr)_+.
\]
There is also a direct version without the binary-depth restriction.
Choose fixed integers \(h_i\ge1\) with \(\sum_{i\in J}1/h_i\le1\).
Hölder on the normalized character average and (4) give
\[
 E(J,k;\gamma)\ll
 X^{\,k\sum_{i\in J}\alpha_i+
       \sum_{i\in J}(k\alpha_i-w/h_i)_++\varepsilon}. \tag{6}
\]
The gap when the reciprocal sum is below one can be assigned to
the constant function one. Interpolating between adjacent integer
orders gives additional bounds.

These statements cover a root group
\(\prod_{i\in J}x_i^{i/g}\) with \(g\mid i\) on the group.
Use \(e_i=i/g\) for the group itself or \(e_i=i\) after raising
the group to power \(g\). An additional fixed exterior power
can be absorbed into the \(e_i\) in the same way.
If the grouped products range over a smaller set, enlarge it to
all products of unit coordinates in their boxes and count tuples.
If other coordinates have already been frozen in each copy,
their monomial products are units and contribute only \(\gamma\).
Freezing them has cost equal to their number of copies times
their total coordinate mass. Thus the estimates apply uniformly
in the boxed twists needed in the collision argument.
[Zeng, Corollary 3.2, Lemma 3.3 and Proposition 4.1](https://hexagonmath.org/pdf/2611.00621v1)

### Uniform energy bounds from interpolation

Write \(y=L_T\), \(w\) for the modulus exponent, and \(\alpha\) for one
coordinate mass. In normalized Hölder take reciprocal weight
\(t=\alpha/y\) and interpolate the logarithm of the norm between
\(1/(m+1)\le t\le1/m\). For fixed \(t\) the penalty beyond the diagonal
mass \(\alpha\) is the chord, as a function of the reciprocal weight,
between \((\alpha-w/(m+1))_+\) and \((\alpha-w/m)_+\). The weights
sum to one. For a fixed dimension, coordinates below an arbitrarily
small mass threshold can be frozen in both copies, with a vanishing
loss and a uniform upper bound on the remaining moment orders.

Put \(q=w/y\). If the chord crosses the kink, its penalty divided by
\(\alpha\) is
\[
 (1-q/((m+1)t))(m+1)(1-mt).
\]
Its maximum is at most \((\sqrt{m+1}-\sqrt{mq})^2\). On the active
range \(m/(m+1)\le q\le(m+1)/m\), this bound decreases with \(m\) at
fixed \(q\). A chord wholly in the positive range gives the relative
penalty \((1-q)_+\), which is also at most
\((\sqrt{2}-\sqrt q)_+^2\). Thus the full residue energy exponent \(K_T\)
obeys
\[
 K_T\le y+(\sqrt{2y}-\sqrt w)_+^2. \tag{7}
\]
The case \(y=0\) follows by the trivial bound on the number of
tuples. At \(w\ge2y\) the exponent is at most \(y\).

The same choice \(t=\alpha/y\) gives a convex upper bound in \(w\).
At \(w=0,y,2y\) it is at most \(2y,(4-2\sqrt2)y,y\), respectively.
Linear interpolation in \(w\), using \(4-2\sqrt2\le47/40\), yields the
rational bound
\[
 K_T\le \max\{2y-\tfrac{33}{40}w,\;
               \tfrac{27}{20}y-\tfrac7{40}w,\;
               y\}. \tag{8}
\]
It uses only the full mass of the released letter and the lower
height of the frozen modulus. One may thus keep the higher
coordinate masses aggregated in the exponent program.

### Using higher moments on both released letters

Fix a residue lock \(O_C\) on \(C\) with mass \(f_C\) and modulus
lower exponent \(w_C\). Choose \(J_T\) and integers \(h_{T,i}\)
as in (6) on each released letter \(T=A,B\).
Freeze both copies of the coordinates outside \(J_T\) when
counting the congruence energy. The bound at \(k=1\) gives
\[
 K_T=2L_T-\sum_{i\in J_T}\min(\alpha_{T,i},w_C/h_{T,i})
\]
as an upper exponent for that energy, up to \(\varepsilon\).
The signed equation makes the residue of one complete value a
fixed unit multiple of the residue of the other.
The same residue Cauchy–Schwarz argument as above gives the
box exponent
\[
 D\le f_C+\frac{K_A+K_B}{2}+\varepsilon. \tag{9}
\]
The analogous substitution into the grouped collision inequality
uses \(w_B\) and the mass of the released coordinates of \(A\)
in place of its full mass. In a projected model the unexposed
higher coordinates may stay outside \(J_T\), at the usual
freezing cost. A fixed projection has only finitely many other
coordinates in the underlying anatomic box.

At the asymmetric profile, work in units of \(1/95920\).
Take \(C=c\) and \(O_C=\{4\}\).
Then \(f_C=11335\) and \(w_C=45340\).
For both free letters take \(J_T=\{1,2,3,4,5\}\) and
\((h_{T,1},\ldots,h_{T,5})=(10,8,4,3,25)\).
Their reciprocal sum is \(509/600\).
In every case \(h_{T,i}\alpha_{T,i}\le45340\), so
\(K_A=31714\) and \(K_B=32103\). Equation (9) gives
\[
 D\le \frac{86487}{191840}+\varepsilon.
\]
One can also lock \(\{3,5\}\) on the short letter \(a\).
The mass is \(10487\) and the modulus exponent \(33993\).
For both free letters choose orders \((8,6,3,3,25)\)
with reciprocal sum \(599/600\).
Only degrees \(1\) and \(4\) incur a penalty, and they
have total penalty
\(4504-33993/8+11335-33993/3=2071/8\) per free letter.
The resulting exponent is \(342791/767360+\varepsilon\).

### Counting two omitted powers with a partial modulus

On \(C\), freeze a coordinate set of mass \(f_C\), whose product
with the fixed coefficient gives a modulus \(M\gg X^{w_C}\).
On each of the other letters \(T=A,B\), select a nonempty group
\(I_T\) whose degrees are all divisible by a fixed \(g\ge1\).
The selected groups are on distinct letters. Their root values are
\[
 u_T=\prod_{i\in I_T}x_{T,i}^{i/g}.
\]
Write \(p_T=\sum_{i\in I_T}\alpha_{T,i}\), and take any upper
exponent \(H_T\) for the size of \(u_T\).
Freeze the complementary coordinates on \(A\) and \(B\).
Writing \(P_A,P_B\) for their products with the respective fixed
coefficients and signs, the additive relation gives
\[
 (u_A/u_B)^g\equiv -P_B/P_A\pmod M.
\]
Primitivity makes the coefficients and both roots units modulo
\(M\). It also makes \(\gcd(u_A,u_B)=1\), since each root divides
the corresponding complete value. The power congruence gives
\(O_{g,\varepsilon}(X^\varepsilon)\) admissible residue classes \(r\)
for \(u_A/u_B\), by the prime-power unit groups and the Chinese
remainder theorem.

For one class \(r\), let \(1\le u_A\le U\), \(1\le u_B\le V\)
with \(U\ll X^{H_A}\) and \(V\ll X^{H_B}\).
Choose a pair \((u_0,v_0)\) in the class with largest numerator.
For every pair \((u,v)\) in it, \(uv_0-u_0v\) is divisible by \(M\)
and has absolute value at most \(UV\). Two pairs with the same
determinant differ by an integer multiple of \((u_0,v_0)\), since
that vector is primitive. Their first coordinates both lie in
\([1,u_0]\), so they coincide. This gives \(O(1+UV/M)\) pairs
uniformly in \(U,V\), including \(M=1\).

The two roots determine the two complete free values, and the
signed equation then determines the third complete value.
The divisor bound controls the number of choices of coordinates
representing each of the three values. Restoring the frozen choices
gives the uniform box bound
\[
 D\le f_C+L_A+L_B-p_A-p_B+(H_A+H_B-w_C)_++\varepsilon. \tag{10}
\]
The two root groups need only have a common fixed divisor
of their degrees. One may conservatively freeze every
unexposed higher coordinate. If the tail of a large \(C\) is
entirely frozen, the height complement can give the same
improved lower bound for \(w_C\) as in the other residue rows.
The reduced-fraction argument is the rectangular form of the
angular reconstruction in the full-letter graph above.
[Powerful-number graph, codegree estimate](POWERFUL_SUM_GRAPHS.md#codegrees-for-any-omitted-exponent)

### A symmetric profile separating the two bounds

Take each letter to have degree masses
\[
 \alpha_1=\frac{19}{150},\qquad
 \alpha_4=\frac{30}{150},\qquad
 \alpha_{11}=\frac1{150},
\]
and zero mass in the other degrees.
Each total mass is \(1/3\) and each weighted height is one.
For the energy methods, consider the ideal substituted exponent
\(K=2s-\min(s,w)\), where \(s\) is the free mass and \(w\) is
the frozen modulus exponent. It is at least as strong as the
moment exponents above. In fact the sum of the savings in
(6), at \(k=1\), is at most \(\sum_i\alpha_i\) and at most
\(w\sum_i1/h_i\le w\). Interpolated moment bounds satisfy the
same two constraints. In a collision row, freezing additional
coordinates on the energy side cannot improve the ideal energy
contribution. The sum of their freezing cost and half of the
ideal exponent is \(L-\tfrac12\min(s,w)\), which is minimized
at the largest free mass.

It suffices to enumerate the eight masks on \(\{1,4,11\}\),
including the degree-\(11\) coordinate individually. Let
\(m(I),h(I)\) be the mass and weighted height of a mask. The
row minima are
\[
\begin{array}{c|rrrrrr}
 \text{row}&\text{pair}&\text{Fourier}&\text{geometry}&
 \text{ideal collision}&\text{ideal direct}&\text{shift}\\
 \hline
 \text{minimum}&2/3&17/30&3/5&8/15&8/15&169/300 .
\end{array}
\]
For the ideal direct bound, the mass and lower height of
\(O_C=\{4\}\) are \(1/5\) and \(4/5\). It attains the
minimum \(8/15\).
For the ideal collision bound, the minimum uses \(I_B=\{4\}\).
Its difference cap is zero and its complementary modulus
has exponent \(1/5\), so \(K=7/15\).
With \(f_B=2/15\) and \(L_C=1/3\), the off-diagonal term is
\(8/15\). The exact scan takes all free masks, all target masks,
and all their subgroup locks and includes the degree-\(11\)
difference cap. The full-letter common-power graph also has
minimum \(8/15\), with just degree \(4\) omitted on both free
letters. Every allowed common-power group removes at most
\(1/5\) of each free letter's mass, and the two degree-\(4\)
roots have combined height \(2/5<1/2\).

The stronger bound (10) instead uses \(O_C=\{4\}\) and
\(I_A=I_B=\{4\}\). Its power is \(g=4\), its two root heights
are \(1/5\) each, and their sum is below \(w_C=4/5\).
Thus (10) gives
\[
 D\le\frac15+2\left(\frac13-\frac15\right)+\varepsilon
     =\frac7{15}+\varepsilon.
\]
This \(8/15\) obstruction concerns the ideal energy and the
other listed rows without the partial-modulus rational count.

### A symmetric obstruction including the partial-modulus rational count

Give all three letters the profile
\[
 (\alpha_1,\alpha_3,\alpha_4)=\frac1{12}(1,1,2),
 \qquad \alpha_i=0\quad(i\notin\{1,3,4\}).
\]
Every letter has mass \(1/3\) and weighted height one.
The point \(D=1/2\) survives both the ideal moment constraints
and (10), along with the earlier geometric, collision, and shift
bounds. The exact row minima are \(2/3\) for the pair bound,
\(7/12\) for the Fourier bound, \(1/2\) for geometry, \(1/2\)
for the ideal direct energy, \(1/2\) for the ideal collision
energy, \(5/8\) for the paired shift, and \(1/2\) for the rational count.

For example, locking degree \(4\) and releasing just degree \(4\)
on both free letters makes the rational row attain \(1/2\).
The degree-one version of the rational row also permits \(1/2\).
In units of \(1/12\), any lock has \(w\le4f\), so
\(f+(H_A+H_B-w)_+\ge (H_A+H_B)/4\). On each letter, the root
mass minus one fourth of its weight is at most one, attained
by releasing degrees \(1,3\). The resulting lower bound on the
right-hand side is \(8-1-1=6\).
[Exact check of the two symmetric profiles](verify_abc_higher_moment_obstructions.py)

### A uniform analytic exponent of \(4/7\)

Put \(d_0=4/7\). In a height shell the largest letter and the
larger addend both have size comparable to the shell height.
Label these two letters \(B,C\), and label the remaining letter
\(A\). The weighted heights of \(B,C\) can be taken to be one,
and the weighted height of \(A\), denoted \(h_A\), is at most
one. Small shell and exceptional-coefficient errors are
absorbed in the final exponent loss. Keep all original
coordinates of the fixed anatomic decomposition.

Write \(L_T=\sum_i\alpha_{T,i}\), \(L=\sum_T L_T\le1\), and set
\[
 x=1-d_0=3/7,\qquad z=2d_0-1=1/7,\qquad F=4-6d_0=4/7.
\]
We may suppose the pair and Fourier bounds do not already
give \(D\le d_0\). Their consequences are
\[
 z\le L_T\le x,\qquad L_B+L_C\ge d_0,\qquad
 \sum_T F_T\le F,\qquad F_T\ge\max_{i\ge2}\alpha_{T,i}.
 \tag{11}
\]
For the first two statements, each mass is at most \(L-d_0\).
Subtracting the other two shows each mass is at least
\(2d_0-L\ge z\). The geometry row for released mass \(s\) and
released weight \(w\) is \(D\le\max(L-s,L-1+w-s)\).

Use \(P(y,w)=(\sqrt{2y}-\sqrt w)_+^2\). The full released energy
on a letter of mass \(y\) is at most \(y+P(y,w)\) by (7),
and is always at most \(2y\). The function \(P\) is increasing
in \(y\), decreasing in \(w\), and jointly convex. For convexity,
\[
 P(y,w)=\sup_{0<t\le1}\{2(1-t)y+(1-1/t)w\}.
\]
All interpolation penalties decrease with the available weight.
For \(q=w/y\ge1\), a chord with index \(m\) can contribute only
when the coordinate mass exceeds \(w/(m+1)\). The relative
penalty on chords with \(m\ge m_0\) is bounded by
\((\sqrt{m_0+1}-\sqrt{m_0q})_+^2\). The absolute penalty of a
single chord with index \(m\) is at most
\(y[m+1-mq]_+^2/(4m(m+1))\), by maximizing the product of
the two linear factors in the chord. Direct residue counting
by (9) gives \(D\le f+(K_U+K_V)/2\) for a lock of mass
\(f\) on the third letter.

Releasing only a single coordinate of degree \(i\ge2\) and
mass \(b\) on \(T\) gives zero difference-fibre penalty in the
grouped collision estimate (1). A nonzero difference of
\(i\)-th powers determines their gap up to divisor-many
choices, and then determines the lower base. Thus, for either
ordering of the other two letters \(U,V\), with \(f=L_T-b\),
\[
 D\le f+\max\{L_V,(L_V+K_U)/2\}. \tag{12}
\]
The modulus in \(K_U\) has weight equal to the frozen
weighted mass on \(T\).

**Large coordinates and early degrees.** Set
\[
 c_0=2/7,\qquad h=.152,\qquad u=.09.
\]
If any degree-\(i\) coordinate with \(i\ge2\) has mass \(b\ge c_0\),
take \(U\) to have the smaller other mass in (12), and use
\(K_U\le2L_U\). The two branches are at most
\((3+x)/4-c_0\) and \(1-z-c_0\), respectively, both equal to
\(d_0\).

If a tall letter has degree-one mass \(b\ge h\), lock all its
other coordinates. The mass term in (9) is at most
\((1+x)/2-b\), and its weight is \(1-b\). Splitting at \(b=.25\)
bounds the exponent by the larger of
\[
 (1+x)/2-h+P(x,.75),\qquad (1+x)/2-.25+P(x,1-x).
\]
This is below \(.56587\).

We can also settle short degree-one mass \(a\ge c_0\) when
\(h_A\ge .5\). Lock its complement at weight at least \(.5-a\).
The two tall masses sum to at most \(1-L_A\), and each is at
most \(x\). By monotonicity and convexity in the masses, the
sum of their penalties is at most the sum at the masses
\(x,1-L_A-x\). Therefore
\[
 D\le (1+L_A)/2-a+
 \tfrac12\{P(x,.5-a)+P(1-L_A-x,.5-a)\}.
\]
This is convex in \(L_A\in[a,x]\). At its two endpoints the
right-hand sides are convex in \(a\in[c_0,x]\), by joint
convexity of \(P\). The resulting endpoint bounds are below
\(.539\).

Now let a tall letter \(T\) have a degree-two or degree-three
coordinate of mass \(b\in[h,c_0]\). Use (12) with \(U\) the
larger other mass. The diagonal term is at most
\((1+x)/2-b\), and the available weight is at least
\(w_0=1-3b\). For \(b\le.175\), all remaining coordinate
masses may be assumed at most \(c_0\). The only case unhandled
so far would be a larger short degree-one mass when
\(h_A<.5\). Then release all of the short letter and the
chosen early coordinate. The two geometry branches are
at most \(1-c_0-h<d_0\) and
\(.5-c_0+2(.175)<d_0\).

For \(b\le.175\), put
\[
 g_m(w)=(\sqrt{m+1}-\sqrt{mw/x})_+^2 .
\]
At weight \(w=w_0\), \(w\ge x\). Without an \(m=1\) penalty
the total released penalty on \(U\) is at most \(xg_2(w)\).
An \(m=1\) penalty requires a coordinate of mass \(a>w/2\).
It is at most \(2(a-w/2)(1-a/x)\), which is increasing for
\(w/2\le a\le c_0\) throughout this range of \(w\).
Its upper bound at \(c_0\) is
\(2(c_0-w/2)_+(1-c_0/x)\). The other total mass is at most
\(x-w/2\) and has no \(m=1\) contribution. For
\(b\le.161\) we have \(w\ge6x/5\), so this other mass is at
most \(w/3\) and contributes only with \(m\ge3\).
Thus the off-diagonal bound is at most
\[
 (1+x)/2-b+
 \tfrac12\max\bigl\{xg_2(w),\,
 2(c_0-w/2)_+(1-c_0/x)+(x-w/2)g_m(w)\bigr\},\qquad w=1-3b,
\]
with \(m=3\) on \([h,.161]\) and \(m=2\) on \([.161,.175]\).
These functions are convex on their intervals. In particular,
the last product is convex because both factors are
nonnegative and decreasing, the first affine and the second
convex. Their endpoint bounds are below \(.56705\).

For \(.175\le b\le 11/42\), use the convex bound
\((1+x)/2-b+P(x,1-3b)/2\), whose endpoint values are below \(.568\).
For \(b\ge11/42\), \(w_0\le x/2\). A released letter of mass
\(y\ge2w_0\) has interpolation penalty \(y-w_0\) at this
weight. Both endpoints of every proportional interpolation
chord are in the linear range. For \(y<2w_0\), (7) gives
penalty at most \(P(2w_0,w_0)=w_0\). In either case the
penalty is at most \(x-w_0\). The off-diagonal exponent is
at most \(x+b/2\le d_0\).

**Two tall letters with large high-degree coordinates.**
Suppose each has a coordinate of degree at least four and
mass \(b_B,b_C\ge h\). Write \(S=b_B+b_C\le.5\) and
\(G(b)=b+P(x,4b)\). A singleton lock on \(T\) gives
\(D\le .5+G(b_T)-L_T/2\). If both bounds exceed \(d_0\), then
\[
 G(b_B)+G(b_C)>2d_0-1+(L_B+L_C)/2\ge 3/7 .
 \tag{13}
\]
Here \(G\) is convex and increasing on \([h,.25]\), since
\(h>8x/25\). Also \(L_T\le1-3 b_T\), because every
remaining coordinate has weight at least its mass.

These conditions force \(S>.42\). If \(S\le.42\) and both
masses are at most \(.228\), their \(G\)-sum is at most
\(G(.192)+G(.228)<.42245\) by convexity and monotonicity.
If one exceeds \(.228\), the other is below \(.192\).
The tall mass sum is then at most
\(1-3(.228)+2G(.192)-z<.56204<d_0\), contradicting (11).

If \(S\ge x\) and at least one of the two degrees is at
least five, release everything except the chosen coordinates.
Their excluded weight minus mass is at least \(3S+h>2-d_0\),
and their excluded mass is at most \(.5\). Both geometry
branches are at most \(d_0\). If \(S\ge x\) and both degrees
are four, use (10) for \(h_A\ge x\), with the whole short
letter locked and the two degree-four coordinates released.
This gives \(D\le\max(L-S,L-h_A)\le d_0\). For \(h_A<x\),
release everything except the chosen coordinates in the
geometry row. The released weight is at most
\(h_A+2-4S<2-3x<1\), so \(D\le S\).

It remains to consider \(.42<S<x\). The singleton failures
and convexity bound \((L_B+L_C)/2\) by
\(G(.25)+G(x-.25)-z<.29223\). Also \(F_A\le F-S<F-.42\).
Put \(\tau=(4d_0-2)/3=2/21\). If a coordinate on \(A\)
of degree at least four has mass \(a\ge\tau\), its singleton
lock gives at most \(.29223+a+P(x,4a)\), a convex expression
on \([\tau,.25]\) whose endpoint bounds are below \(.54224\).
A degree-three mass \(a\ge .08\) likewise gives at most
\(.29223+a+P(x,3a)\), with both endpoint bounds on
\([.08,F-.42]\) below \(.56226\).

Assume therefore that all degree-at-least-four masses on
\(A\) are below \(\tau\) and its degree-three mass is below
\(.08\). Set \(r=2-d_0-3S>0\).
When the short degree-at-least-four coordinates have total
excess, meaning weight minus mass, at least \(r\), omit
them greedily up to the first crossing of \(r\). Their total
mass is at most \(r/3+\tau=d_0-S\). If their excess is below
\(r\) but the short degree-three coordinate supplies the rest,
omit all those higher coordinates and the degree-three
coordinate. Their mass is at most \(r/3+.08<d_0-S\).
With these omissions in addition to the two tall exclusions,
the second geometry branch is at most \(2-3S-r=d_0\).

If the entire short excess from degrees at least three is
below \(r\), the short height satisfies
\(h_A\le x+F-S+r\) because \(\alpha_{A,2}\le F_A\).
Release everything except the two chosen tall coordinates.
The second branch is at most
\(h_A+1-3S\le x+F+3-d_0-7S<.489\), and the first is \(S\).

**A tall letter with only small coordinates.** We are left
with at least one tall letter all of whose masses are below
\(h\). Take \(C\) of greatest mass among those tall letters.
The other tall letter \(B\) has all masses at most \(.25\)
and at most one mass at least \(h\), since all such masses
must now have degree at least four and \(9h>1\).

If the third-largest mass of \(C\) is at least \(u\), let
\(i<j<k\) be the degrees of its three largest coordinates.
For \(k\ge6\), lock that coordinate. Its mass is at most
\(h\), its weight at least \(6u\), and \(L_C\ge3u\).
This gives \(D\le .5+h-3u/2+P(x,6u)<.55348\).

For \(k\le5\), lock all coordinates except degrees \(i,j\),
whose masses we denote by \(a,b\).
The height inequality gives \(ia+jb\le1-ku\) and
\[
 L_C\le M(a,b)=\min\{x,1-(i-1)a-(j-1)b-(k-1)u\}.
\]
Indeed the degree-\(k\) mass is at least \(u\) and all other
weights are at least their masses. The direct-residue
exponent is at most
\[
 .5+M(a,b)/2-a-b+P(x,1-ia-jb).
\]
On each side of the affine boundary in \(M\), this is convex.
Check the vertices of the polygon
\([u,h]^2\cap\{ia+jb\le1-ku\}\) split by that boundary.
All nonempty cases with \(1\le i<j<k\le5\) have bounds
below \(.56822\).

If the third-largest mass is below \(u\), write \(b\) for
the largest and place coordinates in descending mass order
in whichever of two bins is lighter. The final mass
difference is at most \(\max\{u,2b-L_C\}\). Before the first
sign change, the running difference is the largest mass
minus the sum of the others assigned. A sign change first
occurs at the third coordinate or later and all subsequent
differences are at most that coordinate's mass.

If the difference is at most \(u\), choose the bin of
weight at least \(.5\). Its mass is at most \((L_C+u)/2\).
At weight \(.5\) write \(g_2=g_2(.5)\), \(g_3=g_3(.5)\).
A contributing \(m=1\) chord on \(A\) is unique and requires
mass larger than \(.25\). If it is absent, the penalty is
at most \(xg_2\). If its mass is at most \(.28\), its
penalty is at most \(2(.28-.25)(1-.28/x)\) and the rest is
at most \((x-.25)g_2\). If its mass is at least \(.28\),
the penalty is at most \((x-.25)^2/(2x)\) and the other
mass contributes only with \(m\ge3\), at most
\((x-.28)g_3\). These give \(P_A<.03969\).
On \(B\), \(m=1\) is absent and an \(m=2\) penalty occurs
at most once, as it requires mass greater than \(1/6>h\).
The penalty is at most the larger of
\(xg_3\) and \(x(3-1/x)^2/24+(x-1/6)g_3\), giving
\(P_B<.01231\). By (9),
\(D<.5+u/2+(.03969+.01231)/2=.571<d_0\).

If \(2b-L_C>u\), then \(L_C<2h-u=.214<d_0/2\).
Both tall letters cannot have all masses below \(h\) by the
choice of \(C\) and (11). Thus \(B\) has a degree-at-least-four
coordinate with mass at least \(h\) and
\(L_B>d_0-2h+u\). Its singleton lock gives
\(D\le .5+G(.25)-L_B/2=.75-L_B/2<.57129\).

The [constant check](verify_abc_analytic_four_sevenths_constants.py)
verifies the numerical inequalities using exact rational bounds
on the radicals. The case
thresholds and inequalities are stable under an arbitrarily
small perturbation of shell weights, yielding the indicated
\(X^\varepsilon\) loss after shell summation and the fixed
anatomic decomposition. Tiny released coordinates may be
frozen as in (7) to bound the required moment orders
uniformly.

### A determinant bound at the joint \(1/2\) profile

Consider primitive triples whose coordinates on each letter have
scales \(Y,Y,Y^2\) at degrees \(1,3,4\), respectively, where
\(Y=X^{1/12}\). Write the coordinates as \(x_T,y_T,z_T\).
Fixed exceptional coefficients and signs may be attached to the
three products. Their subpower contributions can be included in
the exponent loss. For this profile the box count satisfies
\[
 D\le \frac{249}{500}+\varepsilon. \tag{14}
\]
The proof separates a divisor range and uses two quartic congruence
estimates for the remaining prime factors.

#### Quartic congruence estimates

Let \(p\) be a prime with \(p\ge Y^{2-12\delta}\), where
\(\delta=1/500\), and let \(\gamma\) be a unit modulo \(p^4\).
Uniformly in \(\gamma\), for positive integers \(x\ll Y\) and
\(z\ll Y^2\) that are units modulo \(p\), we have
\[
 \begin{aligned}
 \#\{(x,z):xz^4\equiv\gamma\pmod{p^4}\}
     &\ll_\varepsilon Y^{21/25+\varepsilon},\\
 \#\{(x,z):xz^{-4}\equiv\gamma\pmod{p^4}\}
     &\ll_\varepsilon Y^{3/4+\varepsilon}.
 \end{aligned} \tag{15}
\]
For the first congruence cut the \(z\) range into intervals of
length \(H=Y^h\) with \(h=29/20\), and for the second use \(h=5/3\).
All intervals have length less than \(p\) for large \(Y\).

Take \(R=100\). For the first congruence use the monomials
\(x^a z^b\) with \(0\le a\le R\), \(0\le b\le4a\).
Group them by \(w=b-4a\), so the monomials at each fixed weight are
\(z^w t^a\) with \(t=xz^4\) and consecutive \(a\).
For the second congruence use \(0\le a\le R\),
\(0\le b\le4(R-a)\) and group by \(w=b+4a\), with \(t=xz^{-4}\).
Take a square matrix evaluating the chosen monomials at any
matching number of solutions in one interval.
In a weight class of length \(\ell\), replace the columns by the
evaluations of
\[
 z^w t^{a_0}(t-\gamma)^j,\qquad j=0,\ldots,\ell-1,
\]
where \(a_0\) is the first \(a\) in that class and an integer
representative of \(\gamma\) is used.
These are integral linear combinations of the original monomials
in the same class, with a triangular change of determinant one.
Since \(z\) is a \(p\)-adic unit, column \(j\) is divisible by
\(p^{4j}\), including when \(w\) or the exponent of \(z\) in \(t\)
is negative. Thus the determinant is divisible by \(p^{4S}\), where
\[
 S=\sum_w\binom{\ell_w}{2}
   =4\binom{R+1}{3}+\binom{R+1}{2}
   =\frac{R(R+1)(4R-1)}6 .
\]
Both monomial families have the same multiset of weight-class
lengths. In the second family the length at \(w\) is
\(1+\lfloor w/4\rfloor\) for \(0\le w\le4R\), and in the first
family the lengths appear in reverse order.

Translate \(z\) by the left endpoint of the interval.
Within each fixed \(a\) the exponents \(b\) run from zero to their
maximum, so this translation also preserves the determinant.
Its absolute value is at most a constant times
\(Y^{A+hB}\), where \(A=\sum a\) and \(B=\sum b\) over the chosen
monomials. Put \(A_0=R(R+1)(8R+7)/6\).
For the first family \((A,B)=(A_0,2A_0)\).
For the second family \((A,B)=(S,2A_0)\).
In units of \(\log Y\), the lower bound for the modulus exponent
of the determinant is \((8-48\delta)S\).
At \(R=100\) its excess over \(A+hB\) is \(53833/5\) for the
first family and \(1633574/15\) for the second.
All the evaluation determinants therefore vanish for sufficiently
large \(Y\). The points in the interval lie on a nonzero polynomial
\(F(x,z)\in\mathbb Z[x,z]\) of total degree at most \(5R\).
If there are fewer points than monomials this already bounds their
number by a constant, or the same conclusion follows from the
linear dependence of the evaluation columns.

Factor \(F\) into primitive irreducible polynomials over
\(\mathbb Q\). A factor of total degree less than five for the
first congruence has nonzero reduction modulo \(p\), and its
substitution \(x=\gamma z^{-4}\) is a nonzero Laurent polynomial.
Indeed the kernel of this substitution is the ideal generated by
\(xz^4-\gamma\), of total degree five. Thus the factor gives
\(O_R(1)\) possible residues \(z\bmod p\). The interval has length
less than \(p\), and \(x\ll Y<p\) is fixed by its residue once
\(z\) is chosen. For the second congruence, the same argument
uses the substitution \(x=\gamma z^4\) and factors of degree less
than four. The kernel is generated by \(x-\gamma z^4\), of
total degree four.

For an absolutely irreducible factor of degree \(d\) above the
respective cutoff, apply the coefficient-uniform Bombieri–Pila
bound in the square of side \(O(H)\) enclosing the translated
interval and the \(x\) range. It gives \(O_{R,\varepsilon}
(H^{1/d+\varepsilon})\) points. A rationally irreducible factor
that is geometrically reducible has \(O_R(1)\) rational points,
since these lie in the intersection of distinct conjugate
components and Bézout's bound applies.
Summing over factors and \(O(Y^{2-h})\) intervals gives
\(Y^{2-h+h/5+\varepsilon}=Y^{21/25+\varepsilon}\) for the first
congruence, and \(Y^{2-h+h/4+\varepsilon}=Y^{3/4+\varepsilon}\)
for the second.
[Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

#### Divisors and small cofactors

Suppose a quartic coordinate on any letter has a divisor
\(m\asymp X^\theta\) with
\(\delta\le\theta\le1/6-\delta\).
Use that letter as \(C\) and release only the two degree-four
coordinates on the other letters. Their combined rooted height
is \(1/3\), and fixing the remaining coordinates there costs
\(1/3\). One can either lock \(m^4\) alone or lock \(m^4y_C^3\),
fixing \(m\) and, in the latter case, \(y_C\).
The two resulting rational-lock bounds, before the exponent
loss, are
\[
 \max\left\{\frac13+\theta,\frac23-3\theta\right\},
 \qquad
 \max\left\{\frac5{12}+\theta,\frac12-3\theta\right\}.
\]
Use the second for \(\theta\le1/16\) and the first otherwise.
Both are at most \(1/2-\delta\) on their stated portions.
Choosing one divisor for each tuple and summing over its dyadic
sizes contributes only a logarithmic loss.

In the remaining tuples every quartic coordinate has the form
\(z_T=p_T s_T\) with \(p_T>X^{1/6-\delta}\) prime and
\(s_T\ll X^\delta\). Indeed a prime factor at least \(X^\delta\)
must be larger than \(X^{1/6-\delta}\). If every prime factor is
smaller than \(X^\delta\), a first partial product of those
factors reaching \(X^\delta\) is at most \(X^{2\delta}\), a
divisor in the excluded range. The large prime is unique for
large \(X\) because two such factors would have product larger
than the quartic coordinate.

#### The graph at a fixed quartic coordinate

Fix \(z_C\), and freeze \(s_A,s_B\), which costs \(O(X^{2\delta})\)
choices. Put \(M=z_C^4\asymp Y^8\), and absorb the fixed signs,
cofactors and exceptional coefficients into
\(P=c_A x_A y_A^3\) and \(Q=c_B x_B y_B^3\).
The remaining equation implies
\[
 P\alpha^4+Q\beta^4\equiv0\pmod M,\qquad
 \alpha=p_A,\ \beta=p_B .
\]
There are \(O(Y^2)\) possible \(P\) and likewise \(Q\).
They have fixed signs and heights at most
\(Y^{4+48\delta+o(1)}\), smaller than \(M\).
All these products and primes are units modulo \(M\) on any
edge, and the two primes of a witness are distinct by
primitivity.

Build a simple graph on the coefficients \(P,Q\), using a
witness pair \((\alpha,\beta)\) on each edge. For fixed \(P,Q\)
the congruence determines \(Y^{o(1)}\) residue classes for
\(\alpha/\beta\). The angular determinant bound for coprime
pairs of height \(O(Y^2)\) gives \(O(1)\) pairs in each class.
The two complete values then determine the value on \(C\).
Every such value determines divisor-many tuples. This bounds
the original count by the number of graph edges up to
\(Y^{o(1)}\).

Fix distinct \(P_1,P_2\) and consider their common neighbours
\(Q\), with the corresponding edge witnesses
\((\alpha_i,\beta_i)\) for \(i=1,2\).
Dividing the two congruences shows
\[
 \left(\frac{\alpha_1\beta_2}{\alpha_2\beta_1}\right)^4
 \equiv P_2/P_1\pmod M.
\]
The reduced fraction has height \(O(Y^4)\), so the angular
bound again gives \(O(1)\) fractions per unit root class, since
\(M\asymp Y^8\). Thus there are \(Y^{o(1)}\) fractions in all.

If both prime pairs differ within their letter, the displayed
numerator and denominator are coprime. The cross-letter
primitivity applies to \(\alpha_i,\beta_i\), and the
two within-letter differences give the other two needed
coprimalities. Factoring the two products recovers the four
primes up to bounded multiplicity. Then the first edge
congruence fixes \(Q\bmod M\), and its height and fixed sign
give at most one \(Q\). If both primes repeat within their
letter, the congruences instead give \(P_1\equiv P_2\pmod M\).
Their heights and fixed sign rule this out.

If just the \(A\) prime repeats, the recovered fraction is
\(\beta_2/\beta_1\), which determines both distinct \(B\) primes.
The first edge now gives \(Q\alpha^{-4}\equiv\kappa\pmod{p_C^4}\)
for a fixed unit \(\kappa\). Freeze \(y_B\), with \(O(Y)\) choices.
The second estimate in (15) bounds the remaining choices of
\(x_B,\alpha\) by \(Y^{3/4+\varepsilon}\).
If just the \(B\) prime repeats, recover the two \(A\) primes
instead and get \(Q\beta^4\equiv\kappa\pmod{p_C^4}\).
Freezing \(y_B\) now gives the first estimate in (15), at most
\(Y^{21/25+\varepsilon}\) choices.
The common-neighbour count is therefore
\(O_\varepsilon(Y^{46/25+\varepsilon})\).

Cauchy–Schwarz for the degrees gives
\(O_\varepsilon(Y^{3+(46/25)/2+\varepsilon})
=O_\varepsilon(Y^{98/25+\varepsilon})\) edges.
Restoring the \(O(Y^2)\) choices of the full \(z_C\) and the
\(O(X^{2\delta})\) choices of the two cofactors gives the exponent
\[
 \frac{98/25+2}{12}+2\delta=\frac{373}{750}
\]
on this stratum. The divisor stratum gives the larger
\(249/500\), proving (14). The strict determinant and counting
inequalities also tolerate small perturbations in the box
scales, with a correspondingly small exponent loss. Unselected
coordinates can first be frozen. Their weighted height is small
in a small neighborhood of the profile because the three selected
degrees already account for nearly the full height.
[Exact determinant and exponent checks](verify_abc_quartic_determinant.py)
