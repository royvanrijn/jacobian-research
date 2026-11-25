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

### A determinant bound for symmetric three-coordinate profiles

Fix integers \(2\le j<k\) with \(k\ge4\). On each
letter suppose the coordinates of degrees \(1,j,k\) have sizes
\(Z^\alpha,Z^{1-\alpha},Z\), respectively, with \(0\le\alpha\le1\).
Write \(Z=X^b\), where
\[
 b=\frac1{k+\alpha+j(1-\alpha)}.
\]
Thus the total box mass is \(6b\). With
\(\eta=1/(96k^2)\) the count exponent satisfies
\[
 D\le b(3-\eta)+\varepsilon.
\]
For \(\alpha\ge1/2\) one can take \(\eta=1/56\); for
\(j=3,k=4,\alpha=1/2\) the sharper count is
\(D\le 87/176+\varepsilon\). At the two-coordinate endpoint
\(j=3,k=4,\alpha=0\) one can take \(\eta=1/200\), which gives
\(D\le599/1400+\varepsilon\).
The conclusions allow fixed exceptional coefficients with subpower
heights. All exponents before the final restoration of \(b\) are
in units of \(\log Z\).

#### A prime-power fibre in a rectangle

Let \(m,q\) be coprime positive integers, \(0<a\le1\), and
\(p\ge Z^\rho\), where \(0<\rho\le1\) and \(k\rho>1\).
For units \(\gamma,u,z\), consider
\[
 u^m z^{\pm q}\equiv\gamma\pmod{p^k},\qquad
 u\ll Z^a,\quad z\ll Z.
\]
Take a fixed rational \(s\ge1\) and a block length \(Z^h\),
where \(0<h<\rho\). The following strict inequality is sufficient
for both congruences:
\[
 ma+qs h< k\rho\left(1-\frac1{3s}\right). \tag{A}
\]
Use all monomials \(u^i z^l\) with \(0\le i\le mR\) and
\(0\le l\le qsR\), taking a sufficiently large fixed integer \(R\)
such that \(qsR\) is integral. Group by \(ml\mp qi\). In each
weight class the consecutive monomials have exponent differences
\((m,\pm q)\). With \(t=u^mz^{\pm q}\), the same triangular column
change as for (15), replacing consecutive powers of \(t\) by
powers of \(t-\gamma\), gives \(p^k\) divisibility for each unit
of within-class column order. The number \(S\) of unordered pairs
of columns of equal weight is exactly
\[
 \sum_{\nu=1}^R (m(R-\nu)+1)(q(sR-\nu)+1)
   =\frac{mq}{2}\left(s-\frac13\right)R^3+O_{m,q,s}(R^2).
\]
Indeed each pair differs by \((m\nu,\pm q\nu)\) for one positive
\(\nu\). The determinant is divisible by \(p^{kS}\).
The rectangle is downward closed in \(l\) at each fixed \(i\),
so translating \(z\) to the start of the block preserves the
determinant. The leading exponent bounding its absolute value is
\[
 \frac{mq}{2}(ma s+qh s^2)R^3+O_{m,q,s,a,h}(R^2).
\]
Inequality (A) forces every corresponding determinant to vanish.
Taking a nonzero relation from minors gives an integral auxiliary
polynomial of bounded degree and coefficient height \(Z^{O(1)}\).
This also holds when the evaluation matrix has fewer rows than
columns, by using minors of its rank.

Write \(z=z_0+w\) with \(w\) in the translated block. The modular
binomial in \(u,w\) is
\[
 u^m(z_0+w)^q-\gamma
 \quad\hbox{or}\quad u^m-\gamma(z_0+w)^q .
\]
It is absolutely irreducible modulo \(p\). On the torus this
follows from the primitive character \((m,\pm q)\), which
extends to an integral unimodular change of monomial variables.
The binomial has no coordinate-axis factor in \(u,z\), and
translation preserves irreducibility.

Factor the auxiliary polynomial into primitive irreducible
integral factors. Their coefficient heights are also \(Z^{O(1)}\).
For fixed degrees the largest coefficient of a product of two
real polynomials is bounded below by a constant times the product
of their largest coefficients, by compactness after normalizing
each coefficient norm to one. Thus an integral factor has height
at most a constant times that of the product.
A reduced factor coprime to the binomial has only boundedly many
common residue pairs by Bézout. Each residue fixes \(z\) in the
block, and the full congruence then fixes \(u\) up to at most \(m\)
unit classes modulo \(p^k\), since \(p>m\) for large \(Z\).
The \(u\) range is shorter than \(p^k\). A rationally irreducible
factor that is geometrically reducible also has only boundedly
many rational points, by intersecting its conjugates.

For a remaining absolutely irreducible factor \(F\), its reduction
is divisible by the modular binomial. Write
\[
 T=\max\{ (Z^a)^i(Z^h)^l : u^i w^l \text{ occurs in }F\}.
\]
Weighted degree is additive under multiplication over a field and
cannot increase upon reduction. The top weighted degrees of the
two modular binomials show that
\[
 \log_Z T\ge ma+qh
 \quad\hbox{or}\quad
 \log_Z T\ge \max\{ma,qh\},
\]
respectively. Heath-Brown's affine rectangle estimate bounds the
count on \(F\) by \(O_\varepsilon(Z^{ah/\log_Z T+\varepsilon})\).
This is Theorem 15 with two variables, followed by Bézout with the
coprime auxiliary polynomials supplied there. The logarithm of the
coefficient height costs only a subpower factor. In particular
either sign has \(O_\varepsilon(Z^{h/m+\varepsilon})\) solutions
per block, and \(O_\varepsilon(Z^{1-h+h/m+\varepsilon})\) overall.
[Heath-Brown, *Counting Rational Points on Algebraic Varieties*, Theorem 15](https://ora.ox.ac.uk/objects/uuid%3Ac5ab1264-e082-4664-8c45-ae6d7b6f81a1/files/m73013a42d9d1586d603c5a67baef1e0d)

For \(m=1\) use the two triangles from (15), with \(q\) in place
of \(4\). With a positive power of \(z\) the leading coefficients
of the sums of \(i,l\) are \(q/3,q^2/6\); with a negative
power they are \(q/6,q^2/6\). The within-class pair sum in either
case has leading coefficient \(q/6\). Thus \(2a+qh<k\rho\) or
\(a+qh<k\rho\), respectively, suffice to force the determinants
to vanish. The weighted \(T\) bound just proved gives at most
\(Z^{a/q+\varepsilon}\) points per block, since
\(\log_Z T\ge qh\) for both signs. For \(q=k\) and
\(p\ge Z^{1-d}\), taking
\(h\) arbitrarily close from below to \(1-d-2a/k\) or
\(1-d-a/k\) gives an overall bound
\[
 Z^{d+3a/k+\varepsilon} \tag{B}
\]
covering both signs. The block lengths in this use are positive
for the \(d\) below.

#### Removing intermediate factors and applying the graph

Release the two \(k\)-th power coordinates of the free letters.
Their complement costs \(2\) in units of \(\log Z\). If the root
on the modulus letter has a divisor of size \(Z^t\), the partial
modulus rational count (10) with the divisor raised to the \(k\)-th
power gives
\[
 \max\{2+t,4-(k-1)t\}+\varepsilon .
\]
Consequently intermediate divisors with
\[
 \frac{1+\eta}{k-1}\le t\le1-\eta
\]
already give the desired exponent \(3-\eta+\varepsilon\).
One may choose such a divisor and sum over its dyadic scale.
In the remaining boxes each \(k\)-th power root is \(p_Ts_T\)
with \(p_T>Z^{1-\eta}\) prime and \(s_T\ll Z^\eta\).
Indeed if no prime itself lies above the divisor window, a first
partial product of prime factors reaching the lower threshold
lies within the window, since
\(2(1+\eta)/(k-1)\le1-\eta\). The exceptional prime is unique.
Put \(d=\eta\) for this split.

Freeze the entire modulus root \(z_C\), with
\(M=z_C^k\asymp Z^k\), and the two free small cofactors
\(s_A,s_B\). Their cost is \(1+2d\). The graph of the preceding
quartic proof has \(O(Z)\) vertices on each side, indexed by
\(P=c_Ax_Ay_A^j\) and \(Q=c_Bx_By_B^j\), with the fixed
cofactors included in \(c_A,c_B\). These vertices have size
\(O(Z^{\alpha+j(1-\alpha)+kd+o(1)})<M\), since \(j<k\) and \(kd<1\) in the small-\(\alpha\)
case. In the \(\eta=1/56\) variant, \(\alpha\ge1/2\) instead gives
\(\alpha+j(1-\alpha)\le k/2\). The angular reconstruction of edge and codegree
fractions from residues works since \(M\asymp Z^k\) and the
largest numerator and denominator are \(O(Z^2)\).
With distinct within-letter primes, factorization of the
reconstructed coprime products determines all four primes up
to bounded multiplicity. Two simultaneous repeats give the
same vertex by its size. In a single-repeat case the codegree
is bounded by a fibre modulo \(p_C^k\) of one of the displayed
signs, with either \(x_B\) or \(y_B\) fixed. This is the same
prime recovery as in the quartic argument, with the power \(4\)
replaced by \(k\). Each fixed edge accounts for only subpower
multiplicity. If the common-neighbour exponent is \(c\le1\),
Cauchy–Schwarz and restoration of the frozen roots give an
exponent at most \(5/2+c/2+2d+\varepsilon\).

When \(\alpha\ge28\eta\), freeze \(y_B\) and use (B) with
\(a=\alpha\). The codegree exponent is at most
\[
 1-\alpha+d+3\alpha/k\le1+d-\alpha/4 .
\]
The resulting exponent is at most \(3-\alpha/8+5d/2\),
which is at most \(3-\eta\) for \(d=\eta\).
This also proves the bound with \(\eta=1/56\) throughout
\(\alpha\ge1/2\).

For \(\alpha\le28\eta\), instead freeze \(x_B\). Put \(g=\gcd(j,k)\).
The fibre \(u^jz^{\pm k}\equiv\gamma\bmod p_C^k\) splits into
at most \(g\) congruences \(u^mz^{\pm q}\equiv\gamma_i\bmod p_C^k\),
where \(m=j/g\), \(q=k/g\), and \(\gamma_i^g\equiv\gamma\).
For large \(Z\) we have \(p_C>g\), so the bound on the number
of unit roots follows by reduction modulo \(p_C\) and Hensel lifting.
If \(m=1\), then \(2\le q\le k/2\). Take block length \(Z^{3/4}\) in the
triangle estimates. Both conditions hold since \(a\le1\) and
\(2+3k/8<k(1-\eta)\), and the fibre exponent is at most
\(1/4+a/q\le3/4\). This is below the bound needed below.
If \(m\ge2\), use the rectangular fibre with \(a=1-\alpha\)
and \(s=2k/3\).
Condition (A) holds with \(h\) arbitrarily close from below to
\[
 h_*=\frac{3}{4k^2}-\frac{3d}{2k}.
\]
Indeed \(ma\le k-1\), \(q\le k\), and
\((1-d)(1-1/(3s))\ge1-d-1/(2k)\).
With \(\eta=1/(96k^2)\) and \(d=\eta\) we have \(h_*>68\eta\).
Both cases give a codegree exponent at most \(1+\alpha-h/2\).
The graph gives exponent at most
\(3+\alpha/2-h/4+2\eta\le3-\eta\) on taking \(h\) sufficiently
close to \(h_*\). For \((j,k,\alpha)=(3,4,0)\), use instead
\(\eta=d=1/200\) and \(h\) arbitrarily close to \(h_*=9/200\).
Using the sharper rectangular codegree \(1-h(1-1/j)\),
the exact saving \(h_*(1-1/j)/2-2d=1/200\) gives the stated
\(599/1400\) exponent after the vanishing losses.

At \((j,k,\alpha)=(3,4,1/2)\) take \(\eta=3/88\) and
\(\delta=\eta/3=1/88\). The additional lock that includes the
entire degree-three coordinate covers \(t\) from \(\delta\)
to \(1/2-\eta\): its bound is
\(\max\{5/2+t,3-3t\}+\varepsilon\).
This overlaps the earlier window and extends the excluded
interval to \([\delta,1-\eta]\). In its complement, the large
prime still has cofactor of size \(O(Z^\eta)\). The cofactor
itself is a divisor, so it is smaller than \(Z^\delta\) up to
bounded factors. We can take \(d=\delta\), including in the
prime lower bound. Apply (B) with \(a=1/2,k=4\) to obtain
codegree exponent at most \(7/8+d\). The graph exponent is
at most \(47/16+5d/2=3-\eta\).
Because \(b=1/6\) in this profile, the count exponent is
\(87/176+\varepsilon\).

#### Different root scales and degrees on the three letters

The same estimate applies when the scale of the common degree-\(k\)
coordinate on letter \(T\) is \(X^{b_T}\), the other degrees are
\(1,j_T\) with \(2\le j_T<k\), and their two masses sum to \(b_T\).
They may split that mass independently on each letter. Put
\(b=\max_T b_T\), choose a rational \(0<\sigma\le1/2\) with
\(j_T\le(1-\sigma)k\) on each letter, and set \(\eta=\sigma^2/96\).
If
\[
 b_T\ge \max\left\{\frac4{k+1}+\eta,\frac{3\sigma^2}4\right\}b
 \qquad\text{on each letter},
\]
then
\[
 D\le b_A+b_B+b_C-\eta b+\varepsilon.
\]
The complete heights may vary. The choice \(\sigma=1/k\)
gives the bound at all lower degrees \(j_T< k\). With
\(b_T=b\) and fixed \(\sigma\), the saving is proportional
to the total mass \(6b\) uniformly in \(k\).

On a modulus letter \(T\), a divisor of top-root mass \(t\) gives
the rational-lock exponent
\[
 \max\{b_U+b_V+t,\ 2(b_U+b_V)-(k-1)t\}+\varepsilon,
 \qquad \{U,V\}=\{A,B,C\}\setminus\{T\}.
\]
Use the divisor window from
\[
 \max\left\{\eta b,\frac{b_U+b_V-b_T+\eta b}{k-1}\right\}
 \quad\text{to}\quad b_T-\eta b.
\]
The displayed lower endpoint is positive and twice it is at most
the upper endpoint, using \(b_U+b_V\le2b\) and
\((k+1)(b_T-\eta b)\ge4b\). The greedy divisor argument therefore
gives a unique prime with cofactor \(O(X^{\eta b})\) on each
remaining top root.

Take \(b_C=b\) and \(b_A\le b_B\). Freeze \(z_C\) and the
two small cofactors. The modulus \(X^{kb}\) still permits the
angular reconstruction, since the largest wedge numerator and
denominator are \(O(X^{b_A+b_B})\). All graph vertices have absolute
height less than the modulus by \(j_T\le(1-\sigma)k\) and
\(\eta<\sigma\).
In units of \(\log Z=b\log X\), define \(f_T=b_T/b\), and let \(a\)
be the degree-one mass on \(B\) in those units. Its other low
coordinate has mass \(f_B-a\). A block of length \(Z^h\) now has
overhead \((f_B-h)_+\), since each free prime is \(O(Z^{f_B})\).
For \(a\ge28\eta\) use the triangles and freeze the degree-\(j_B\)
coordinate. If \(f_B\) exceeds the block length exponent, using
\(f_B\le1\) recovers the codegree exponent
\(f_B+\eta-a(1-3/k)\). If there is only one block, the codegree
is at most \(f_B-a+a/k\), which is smaller.
For \(a\le28\eta\) use the rectangles after freezing the degree-one
coordinate and extracting the common roots as before. With
\(m\ge2\), use \(s=2/(3\sigma)\) in (A). Since
\(m(f_B-a)\le(1-\sigma)k\) and \(q\le k\), we may take \(h\)
arbitrarily close from below to
\[
 h_*=\frac{3\sigma^2}4-\frac{3\sigma\eta}2>68\eta.
\]
We have \(f_B>h_*\), so the codegree exponent is
at most \(f_B+a-h/2\), with \(h\) close to \(h_*\).
If \(m=1\), the triangle with length \(Z^{3/4}\) bounds the
fibre exponent by \((f_B-3/4)_++f_B/2\), which is also at most
\(f_B-h/2\).

The vertex counts are now \(O(Z^{f_A})\) and \(O(Z^{f_B})\).
With codegree exponent \(c_B\), Cauchy–Schwarz gives an edge
bound with exponent
\(\max\{f_B,f_A+(f_B+c_B)/2\}\).
The two codegree bounds give \(c_B\le f_B-6\eta\) and
\(f_A>3\eta\). Restoring the three frozen factors costs
\(1+2\eta\), leading to the stated exponent.

### Separating repeated residues from the diagonal

Fix an anatomic box and lock coordinates of total mass \(f_C\) on
letter \(C\). Let \(M\gg X^{w_C}\) be a divisor of its complete
value determined by the frozen data. Let \(r_T(t)\) count distinct
permitted complete values on a free letter \(T=A,B\) in the
unit class \(t\bmod M\). We align the two residue tables by
the fixed unit twist from the signs and any fixed coefficients.
Each pair of congruent values determines the third complete value,
and that value has divisor-many coordinate representations.
Thus, before restoring the frozen choices, it suffices to count
\[
 N=\sum_t r_A(t)r_B(t).
\]
Write \(E_T^\ne=\sum_t r_T(t)(r_T(t)-1)\) and
\(E_T=\sum_t r_T(t)^2\). Split off the first representative of
each nonempty \(A\) residue class. Cauchy–Schwarz and
\((r-1)^2\le r(r-1)\) for integers \(r\ge1\) give
\[
 N\le \sum_t r_B(t)+\sqrt{E_A^\ne E_B}.
\]
If \(E_A^\ne\ll X^{U_A+\varepsilon}\) and
\(E_B\ll X^{K_B+\varepsilon}\), this gives
\[
 D\le f_C+\max\{L_B,(U_A+K_B)/2\}+\varepsilon .
 \tag{SR}
\]
The labels \(A,B\) may be interchanged.

For either free letter \(T\), select a nonempty group of
coordinates all of whose degrees have a fixed common divisor
\(g_T\ge2\). Write \(p_T\) for its mass and \(h_T\) for an
upper exponent of the complete value. The off-diagonal energy
obeys
\[
 U_T=2(L_T-p_T)+(h_T-w_C)_+ .
\]
Indeed, freeze both copies of the complementary coordinates.
Their factors form fixed coefficients \(P,Q\). The grouped
roots \(x,y\) then satisfy
\[
 Px^{g_T}-Qy^{g_T}=jM,\qquad j\ne0 .
\]
There are \(O(1+X^{h_T-w_C})\) possible values of \(j\).
For each \(j\) the uniform binary-form bound gives
\(O_\varepsilon(X^\varepsilon)\) pairs \((x,y)\).
This uses the nonzero weighted power-difference count,
uniform without coprimality assumptions
[Zeng, Lemma 3.1](https://hexagonmath.org/pdf/2611.00621v1).
For counting distinct complete values there is no additional
multiplicity in recovering them from the two roots.
The congruence count works at any divisor \(M\),
including a partial-coordinate lock chosen in a divisor stratum.

One may use any valid total-energy bound for \(K_B\) in (SR).
In particular, the diagonal contributes at most \(X^{L_B}\),
so \(K_B\le\max(L_B,U_B)\) from the same off-diagonal argument
on \(B\). In the ordered-tall projection use \(h_A,h_B\le1\)
and the available lower weight \(w_C\le1\) of the locked
modulus. Equation (SR) then becomes a maximum of three affine
bounds,
\[
 D\le f_C+\max\left\{
 L_B,\quad
 L_A-p_A+\frac{L_B+1-w_C}{2},\quad
 L_A+L_B-p_A-p_B+1-w_C
 \right\}+\varepsilon .
 \tag{SR-affine}
\]

At the asymmetric unit-budget profile with units \(1/1424\),
the nonzero coordinates are
\[
 (\alpha_{A,1},\alpha_{A,3})=(266,386)/1424,\qquad
 (\alpha_{B,1},\alpha_{B,4})=(200,306)/1424,\qquad
 (\alpha_{C,1},\alpha_{C,7})=(73,193)/1424.
\]
Lock the entire \(C\)-tail, which has exact weight \(1351/1424\)
from the height of \(C\). Use the singleton degrees \(3,4\)
on \(A,B\). Thus \(U_A=605/1424\), \(U_B=473/1424\), and
\(K_B\le506/1424\). Equation (SR) gives
\[
 D\le \frac{193}{1424}+\frac{605+506}{2848}
 =\frac{1497}{2848}<\frac{27}{50}.
\]
The arithmetic inequality is strict on an open neighborhood of
the profile.

### Repeated residues in a short-addend strip

Suppose the locked letter \(C\) is the short addend, so that
the complete free letters \(A,B\) satisfy \(|A-B|=C\).
For the fixed lock data write \(C=M c\), where \(M\) is the
locked modulus, and suppose \(c\ll X^{s_C}\).
One may also freeze coordinates on the free sides, at masses
\(f_A,f_B\). For \(T=A,B\), let \(P_T\) be its fixed factor,
assume \(P_T\gg X^{v_T}\), and write \(y_T=L_T-f_T\).
Select a nonempty common-power group on each free side,
with fixed root power \(g_T\ge2\) and mass \(p_T\).

Partition the complete values of the two free sides into
half-open intervals of length \(H\asymp M X^{s_C}\), chosen
large enough that \(C\le H\). Only intervals with indices
differing by at most one can supply a valid pair. For one of
the three offsets fix temporarily the corresponding pair of
intervals. Each congruence has the same coefficient twist
after cancelling \(P_A,P_B\) modulo \(M\). Apply the
separated-diagonal estimate to the two residue tables in
these intervals. Sum the diagonal term over the intervals,
and sum their Cauchy–Schwarz terms by Cauchy–Schwarz once
more. For a fixed offset each interval occurs exactly once
on its respective side, so the resulting repeat energies are
sums over disjoint intervals and incur no cost for their
number.

Two distinct \(T\)-values in the same interval and residue
class have free products differing by \(jM\) with \(j\ne0\).
This follows since \(P_T\) is a unit modulo \(M\).
The size of the difference of the complete values gives
\(|j|\ll X^{s_C-v_T}\). Freeze both copies of the remaining
free coordinates outside the common-power group, at cost
\(2(y_T-p_T)\). The weighted power-difference estimate gives
only \(X^\varepsilon\) choices for each nonzero shift. Thus
the sum of the restricted repeat energies has exponent
\[
 U_T^{\rm strip}=2(y_T-p_T)+(s_C-v_T)_+ .
\]
If \(s_C<v_T\), these repeat energies are zero for large \(X\).
In general the total energy summed over the intervals on side
\(B\) is at most \(X^{\max(y_B,U_B^{\rm strip})+\varepsilon}\).
The resulting exponent including all frozen choices is
\[
 D\le f_A+f_B+f_C+
 \max\left\{y_B,\frac{U_A^{\rm strip}+y_B}{2},
              \frac{U_A^{\rm strip}+U_B^{\rm strip}}2\right\}
 +\varepsilon .
\]
The two tall letters may be interchanged.

For the simple case without extra freezing, take
\(f_A=f_B=v_A=v_B=0\). Then \(s_C\) is any upper weight of the
complementary coordinates on the short letter. If the lock
contains the entire tail, this is the exact early complement
weight. Otherwise use the upper bound from the full weighted
height at most one minus the known early weight in the lock.
Either choice is an affine overestimate in the six-index
program. This improves the global shift \(1-w_C\) when a
large part of the bound one was unoccupied height of the
short letter.

### Affine curves after locking a frozen cofactor

On each of two letters \(A,B\), release a nonempty subset of
coordinates whose degrees are multiples of a fixed integer
\(p\) and \(q\), respectively. Write the two released products
as \(x^p,y^q\), with \(p,q\ge1\). Write \(m_A,m_B\) for their
total radical masses and \(\beta,\gamma\) for the numerical
height exponents of \(x,y\). Freeze all remaining coordinates,
including the whole third letter. The signed equation is
\[
 Ux^p+Vy^q=W
\]
for fixed nonzero \(U,V,W\), with \(\gcd(U,VW)=1\).
Let \(w_A\) be a lower bound on the numerical weight of
\(\lvert U\rvert\). Then
\[
 D\le L-m_A-m_B+\frac{(\gamma-w_A)_+}{p}.
 \tag{PC}
\]
The symmetric bound uses \((\beta-w_B)_+/q\).

Reduction modulo \(\lvert U\rvert\) confines \(y\) to
\(X^\varepsilon\) residue classes. Indeed a unit has at most
\(q\gcd(q,2)\) roots modulo any prime power, by the structure
of the unit group, and \(q\) is fixed. For each class use a
representative \(y_0\) between \(0\) and \(\lvert U\rvert-1\) and
write \(y=y_0+\lvert U\rvert z\). There are \(O(1)\) possible
\(z\) if \(\gamma\le w_A\), up to the dyadic scale constants.
Otherwise the box for \(x,z\) has exponents
\(\beta,\gamma-w_A\).
The plane curve
\[
 Ux^p+V(y_0+\lvert U\rvert z)^q-W=0
\]
is absolutely irreducible. For \(p=1\) it is linear in \(x\)
with a nonzero constant leading coefficient. For \(p>1\),
apply Eisenstein in \(x\) at any simple zero of
\(V(y_0+\lvert U\rvert z)^q-W\) over \(\overline{\mathbb Q}\).
All of those zeros are simple because \(VW\ne0\).
The weighted affine plane-curve estimate in equation (3) of
Browning–Heath-Brown uses exponent
\(\log P_2\log P_3/\log T\) in a box of sides \(P_2,P_3\).
The monomial \(x^p\) alone gives \(\log T\ge p\log P_2\).
The estimate gives \(X^{(\gamma-w_A)/p+\varepsilon}\)
points per residue class when \(\gamma>w_A\), uniformly in the
fixed coefficients. Vanishing side exponents follow by
enlargement through \(X^\varepsilon\); if \(\gamma\le w_A\),
each of the \(O(1)\) possible \(z\) determines at most one
positive \(x\). A grouped root has divisor-many
factorizations. Summing over all frozen choices proves (PC).
[Browning–Heath-Brown, equation (3)](https://arxiv.org/pdf/math/0506497)

This gives a maximum of two affine bounds when \(w_A\) is
supplied by the ordered-height model. For a tall \(A\) use
\(w_A=1-p\beta\) in the limiting box. For a short \(A\) use the
sum of the known weights on the frozen indices, with
\(7\) times the frozen tail mass as its tail contribution.

At the cubic–quartic witness of denominator \(6274\), take
\(p=3,q=4\), \(w_A=1189/6274\), \(\gamma=1342/6274\), and
\(L-m_A-m_B=3237/6274\). The new exponent is at most
\(1644/3137<53/100\) at this witness.

#### Shortening the plane box by the short height

In (PC) suppose the square or higher power \(x^p\) is
released on the short addend, and the released \(y^q\)
lies on a tall term. Freeze the other tall term as before.
Then the same count gives the stronger bound
\[
 D\le L-m_A-m_B+\frac{(\gamma+p\beta-1)_+}{p}.
 \tag{PCS}
\]
Indeed let \(h\) be the actual short height on the
logarithmic scale. With the fixed coefficients of (PC),
the two tall terms have size comparable to \(X\).
Their difference is the short term, so the mean value
theorem confines \(y\) to an interval of
length \(O(X^{\gamma+h-1})\) about the positive
real root of the equality between the two tall terms.
Consecutive positive values of \(y\) in one residue
class modulo \(\lvert U\rvert\) are separated by
\(\lvert U\rvert\asymp X^{h-p\beta}\).
Choosing the representative closest to the centre
shortens the box for the new integer coordinate to
exponent \((\gamma+p\beta-1)_+\). The shifted curve
is still absolutely irreducible, and its monomial
\(x^p\) still supplies the same lower bound for the
leading scale in the plane estimate. This proves (PCS).
It uses the short cofactor even when the weighted
height carried by its tail is unknown.

## Collinear congruence lattices for repeated power residues

Fix a locked letter factor \(M\) of root mass \(f\) and a
lower numerical weight \(w\), so \(M\gg X^w\). On each
other letter \(T\), choose a group of degree-\(p_T\)
powers, with logarithmic root size \(r_T\), root mass
\(m_T\), and \(p_T\ge2\). Its free complete values have
the form \(P u^{p_T}\), where \(P\) is the product of
their complementary coordinates. The copies of \(P\)
are not fixed in the full application. Each has at most
\(X^{L_T-m_T+\varepsilon}\) choices. Suppose
\[
 2r_T\le w
\]
for each chosen group. Let \(s\ge0\) be an upper exponent
for the nonzero quotients \(j\) in the repeated residues
of either letter, so each such pair satisfies
\[
 P u^{p_T}-P'(u')^{p_T}=jM,\qquad 0<|j|\ll X^s.
\]
Ordinarily one can take \(s=1-w\) by the size of complete
letters. When \(M\) is a lock from a short addend, pair
the tall letters in intervals of length near
\(M X^{s_{\mathrm{short}}}\), as in the existing strip
argument, and take \(s=s_{\mathrm{short}}\), the upper
weight of the complementary short factors. For a whole
short-tail lock this is the exact early short height.
One may choose the smaller of the two bounds.

Then for every fixed \(M\) the sum of the off-diagonal
energies on \(T\), over the residue classes and any
disjoint intervals used, is at most
\[
 X^{E_T+\varepsilon},\qquad
 E_T=2(L_T-m_T)+s/p_T.
\]
To prove this, fix \(P,P'\). All quantities involved are
units modulo \(M\) when there are admissible points, by
primitivity of each additive triple. The ratio \(u/u'\)
is a \(p_T\)-th root of a unit modulo \(M\) and has
\(X^\varepsilon\) possible residues. This follows by CRT
from the prime-power unit groups. Each possible ratio
\(\lambda\) puts \((u,u')\) in the lattice
\(u\equiv\lambda u'\pmod M\) of index \(M\) in
\(\mathbb Z^2\).

In the dyadic square with root bounds \(U\ll X^{r_T}\),
there are only \(O(1)\) lattice points or all admissible
points are collinear with the origin. Indeed if the box
contains two independent lattice vectors, at least one
of them is independent of a shortest nonzero vector \(v_*\).
The determinant of that pair is a nonzero multiple of
\(M\), so \(\lVert v_*\rVert\gg M/U\gg U\).
Packing points with this separation in the square gives
\(O(1)\). With a strict exponent gap \(w>2r_T\),
two independent vectors cannot occur at all. A vanishing
dyadic exponent slack produces only a freely small power
loss, by the same packing comparison.

In the collinear case write each positive pair as
\((u,u')=d(v,v')\) for one fixed primitive integer ray.
For off-diagonal pairs \(P v^{p_T}-P'(v')^{p_T}\ne0\).
Primitivity of either copy against \(M\) implies
\(\gcd(d,M)=1\). Since
\(d^{p_T}(P v^{p_T}-P'(v')^{p_T})=jM\), it follows
that \(d^{p_T}\mid j\). In particular
\(d\ll X^{s/p_T}\). Count the rays by the number of
congruence lattices, and count \(P,P'\) by the two frozen
complement costs, to obtain \(E_T\). If \(P,P'\) can be
represented by several tuples, the standard divisor
bound absorbs those multiplicities.

Designate one free side as \(B\) and the other as \(A\).
The separated-diagonal inequality gives
\[
 D\le f+\max\left\{L_B,\frac{E_A+L_B}{2},
                        \frac{E_A+E_B}{2}\right\}.
\]
For each residue where \(A\) occupies it once or less, its
pair contribution is at most that of \(B\). In the
remaining residues a multiplicity \(n_A\) satisfies
\(n_A\le\sqrt{2n_A(n_A-1)}\), so Cauchy–Schwarz uses the
repeat energy of \(A\) and the total energy of \(B\).
In the strip version the same reasoning applies for
each of three interval offsets. A second Cauchy–
Schwarz sums over the disjoint intervals in each offset,
without a cost for the number of intervals. Both choices
of \(B\) are valid. The formula also admits the existing
repeat bound on one side where the root-box guard fails.

## A bound without root-box guards

The root-box hypothesis can be removed from the energy formula.
In the rank-two case, choose a shortest lattice vector of
length \(\ell\) and complete it to a basis with second vector's
parallel component of absolute value at most \(\ell/2\).
The perpendicular component has length \(H=M/\ell\), and
\(H\ge\sqrt3\ell/2\) by minimality of \(\ell\).
If the box of side \(O(U)\) contains two independent
vectors, then \(\ell\ll U\), and at least one point has
nonzero perpendicular component, proving \(H\ll U\).
The integer coefficients in this basis are bounded by
\(O(U/\ell)\) and \(O(U/H)\), respectively. This gives
\(O(U^2/M)\) lattice points.

Combining the rank-two count with the collinear count
always gives
\[
 E_T=2(L_T-m_T)+\max\{2r_T-w,\ s/p_T\}.
\]
One can optionally take the minimum of the last maximum
and the old gap \(s\), independently on either side.
In the ordinary residue version \(s=1-w\), the replacement
is already bounded by \(s\), since \(2r_T\le 1\).
With \(p_T\ge3\) and \(s>0\), this ordinary-step improvement
is strict. The bound is continuous and piecewise affine
in the retained scales and needs no auxiliary guard.
Writing the last displayed \(D\) bound as one maximum
of affine functions uses at most seven terms, one for
\(L_B\), two for \((E_A+L_B)/2\), and four for
\((E_A+E_B)/2\). Both displayed numerical examples still
use the collinear value \(s/p_T\) on both sides.


### Exact examples

For the profile with \(D=73/136\), lock the entire short
\(a\)-tail. It has mass \(f=146/1088\) and numerical
weight at least \(7f=1022/1088\), strictly greater than
twice the root exponents of the \(b_4,c_3\) blocks.
The strip quotient weight is \(s=a_1=66/1088\).
Use \(p_B=4,m_B=r_B=212/1088,L_B=372/1088\) on \(b\)
and \(p_A=3,m_A=r_A=292/1088,L_A=504/1088\) on \(c\).
Then
\[
 E_B=\frac{673}{2176},\qquad E_A=\frac{446}{1088},
 \qquad
 D\le\frac{146+409}{1088}=\frac{555}{1088}<\frac{107}{200}.
\]
All comparisons are strict at this profile. The guard
uses only the minimum height of the locked tail, so it
does not require a coordinate for degree seven.

For the profile with \(D=414/781\), lock the \(c_3\)
coordinate of mass \(f=207/781\) and numerical weight
\(w=621/781\). In the ordinary, unrestricted residue
version take \(s=1-w=160/781\) and take the groups
\(a_4,b_5\) of degrees four and five. Both guards are
strict, since \(2r_A=320/781\) and
\(2r_B=207/781\). The complementary masses are
\(205/1562\) and \(113/1562\), respectively. Designate
\(B=b\), so \(L_B=160/781\). The repeat exponents are
\(E_A=245/781\) and \(E_B=145/781\). Thus
\[
 D\le\frac{207}{781}+
    \max\left\{\frac{160}{781},\frac{405}{1562},
                                    \frac{390}{1562}\right\}
   =\frac{819}{1562}<\frac{53}{100}.
\]

### Residue lattices intersected with the real short strip

When the locked term is the short addend, the rank-two
part of the preceding energy bound can use the short
interval itself. Write \(M\) for the fixed modulus
and \(s\) for an upper weight of the complementary
short factors. For repeated residues on either tall
letter the allowable quotient \(j\) is nonzero and
\(\lvert j\rvert\ll X^s\). Under the same notation as
above the new energy bound is
\[
 E_T\le 2(L_T-m_T)+\max(2r_T+s-1,s/p_T).
 \tag{AS}
\]
The rank-two and rank-one cases use the same
congruence classes for \(u/u'\).

To prove this, fix \(P,P'\) and one such class.
Both positive complete terms \(Pu^{p_T},P'(u')^{p_T}\)
are of size comparable to \(X\); both roots are of
size comparable to \(U=X^{r_T}\).
The relation
\[
 \lvert Pu^{p_T}-P'(u')^{p_T}\rvert\ll M X^s
\]
and the mean value theorem give
\(\lvert u-(P'/P)^{1/p_T}u'\rvert\ll U M X^{s-1}\).
Together with the root box this confines the pair to
a strip of area \(O(U^2 M X^{s-1})\).
Apply the rank-two packing calculation after linearly
rescaling this strip to a fixed square. The new
covolume is comparable to the original covolume \(M\)
divided by the area of the strip. If there are two
independent admissible vectors, the resulting bound
is \(O(U^2 X^{s-1})\). If instead the vectors are
collinear, the primitive-ray argument still bounds
their scale by \(O(X^{s/p_T})\). Summing over
coefficient choices and congruence classes gives (AS).
The cancellation in the area-to-covolume ratio uses
the actual modulus size and requires no estimate for
the hidden height of its factors.

For a comparison at unit budget, take short masses
\(a_1=4/25,a_3=1/25\), tall masses
\(b_1=3/25,b_2=1/5,b_4=3/25\), and
\(c_2=3/10,c_6=1/50\) with \(c\)-tail mass \(1/25\)
at limiting weight \(7/25\).
Lock \(a_3\) together with the zero-mass short tail, and release
just the degree-two coordinates on the tall terms, designating \(B=c\).
Here \(s=4/25\) and the lock has weight \(3/25\).
The progression-based strip energies are \(16/25\) and \(7/25\)
on \(b,c\), giving \(D\le27/50\) from the separated
residue formula. Intersecting the lattices with the
strip replaces those energies by \(14/25\) and \(1/5\),
so the same formula gives \(D\le1/2\).
This comparison uses a feasible mass and height
profile for the two energy estimates.

### Asymmetric power groups in the repeat energies

The two copies in an off-diagonal repeat may release different
power groups. Fix a lock of numerical weight at least \(w\)
and mass \(f\), with modulus \(M\). Let \(s\ge0\) bound the
logarithmic magnitude of the nonzero quotient in a repeat
difference divisible by \(M\). The ordinary choice is \(s=1-w\).
The short-addend strip also permits the projected numerical
weight of the unlocked short coordinates.

For a free letter \(T\), choose groups \(I,J\) in the first
and second copies of its repeat. They can intersect or coincide.
Require a common integer \(p\ge2\) dividing every degree in
\(I\cup J\). The groups have factor masses \(m_I,m_J\) and
numerical root exponents \(r_I,r_J\), where the root products
are formed after division by \(p\). Then the repeat exponent
can be taken to be
\[
 E_T(I,J)=2L_T-m_I-m_J+
            \max\{r_I+r_J-w,\ s/p\}. \tag{A}
\]
Choose \(p=\gcd(I\cup J)\) for the strongest form of this
estimate. This form allows larger numerical root exponents
on one copy if the other keeps them below the rank-two
threshold. In the real short strip, the same argument as above
uses a strip with area \(O(UV M X^{s-1})\), so the area term in
(A) becomes \(r_I+r_J+s-1\).

Freeze the complementary products in each of the two
copies, at combined cost \(2L_T-m_I-m_J\). Their complete
values are now \(P u^p,P'(u')^p\). On an admissible nonzero
repeat the difference is \(jM\), with \(1\le |j|\ll X^s\).
All four quantities are units modulo \(M\), since each
complete value is coprime to the locked factor in its
additive triple. Thus \(uu'^{-1}\) lies in one of
\(X^\varepsilon\) classes modulo \(M\), from the CRT bounds
for unit \(p\)-th roots. Each class gives a lattice of
index \(M\) in \(\mathbb Z^2\).

Scale the root rectangle with sides \(U\ll X^{r_I}\),
\(V\ll X^{r_J}\) to a fixed square. In each image lattice,
if the admissible points contain independent vectors, one
can use a Gauss basis with shortest length \(\ell\) and
perpendicular height \(H\). Independence in a fixed square
gives \(\ell,H\ll1\) as in the symmetric proof. The integer
coordinates in this basis obey \(O(1/\ell)\) and \(O(1/H)\)
bounds, so there are \(O((\ell H)^{-1})=O(UV/M)\) points.
If the points are collinear through the origin, each
positive integer pair is a positive multiple \(d(v,v')\)
of the one primitive integer ray. Then \(d^p\mid j\) by
the difference identity and \(\gcd(d,M)=1\). Thus a
collinear residue contributes \(O(X^{s/p})\) points.
Summing over classes and frozen data proves (A). The same
argument aggregates over disjoint strip intervals since
every ordered repeat belongs to at most one interval.

Insert the two \(E\)-bounds in the usual separated
diagonal estimate
\[
 D\le f+\max\left\{
      L_B,\ (E_A+L_B)/2,\ (E_A+E_B)/2\right\}.
\]
The resulting bound is again a maximum of at most seven
affine functions for fixed groups on each side. Legacy
upper bounds on either repeat exponent can also be reused.

#### Exact comparison at a unit-budget profile

Take the following exact nonzero masses:
\[
 a_1=23/600,\quad a_2=3/10,\quad a_6=1/60,\qquad
 b_3=1/25,\quad b_4=29/200,\quad b_5=3/50,\qquad
 c_2=73/200,\quad c_7=1/100,\quad c_8=1/40.
\]
Then \(L_a=71/200,L_b=49/200,L_c=2/5\), the budget is
one, and both \(b,c\) have height one. The last two
coordinates can be combined into the degree-at-least-seven
mass for a projected computation. Lock \(b_5\), so
\(f=3/50,w=3/10,s=7/10\).

On \(a\), choose \(p=2,I=\{2\},J=\{2,6\}\). The two
root exponents are \(3/10,7/20\), and the gap in (A) is
\(7/20\). The resulting \(E_a=133/300\). The best
symmetric choice with this lock has \(E_a=23/50\):
the respective choices of the degree-2 block alone,
the degree-6 block alone, and their union give
\(23/50,\ 119/150,\ 143/300\), with the singleton
degree-6 block counted at \(p=6\). Including the
degree-1 block cannot improve this symmetric
comparison, since (A) requires a genuine power group.

On \(c\), use \(p=2\) and the degree-2 coordinate in
both copies, giving \(E_c=1/2\). This remains the best
symmetric group even when \(c_7,c_8\) are separately
available: the group containing both \(c_2,c_8\)
has rank-two gap \(63/100\) and repeat exponent
\(13/20\), and any group omitting \(c_2\) has
complementary expense at least \(73/100\).
This yields the mixed cap
\[
 D\le319/600<107/200
\]
compared to \(27/50\) for the best symmetric cap with
this same lock. The comparison uses the same fixed modulus throughout.

### Different powers on the two copies

More generally, allow the first and second released groups
to use degrees \(p,q\), respectively. Let
\(g=\gcd(p,q)\), \(P=p/g\), \(Q=q/g\).
The root exponents before lifting are \(r_I,r_J\).
Lift the roots \(u,u'\) to \(Y=u^P,Z=(u')^Q\).
Their root boxes have exponents \(Pr_I,Qr_J\). The same
argument bounds the repeat exponent by
\[
 E_T(I,J)=2L_T-m_I-m_J+
       \max\{Pr_I+Qr_J-w,\ s/\operatorname{lcm}(p,q)\}. \tag{B}
\]
The lock and the nonzero shift satisfy the hypotheses
of (A). Degrees \(p=1\) or \(q=1\) are also allowed,
with the corresponding root just the entire numeric
group product. This allows an aggregate tail in such
a group using its projected upper numeric weight.

The ratio congruence is on \(Y/Z\), with \(g\)-th power
fixed modulo \(M\), including an ordinary single linear
ratio when \(g=1\). The rank-two case therefore gives
\(O(X^{Pr_I+Qr_J-w})\) as before. In a real short strip the rank-two term becomes
\(Pr_I+Qr_J+s-1\). In the collinear case
\((Y,Z)=d(Y_0,Z_0)\) on a fixed primitive integer ray.
As before, \(d^g\mid j\). Since \(Y,Z\) are respectively
\(P\)-th and \(Q\)-th powers, at each prime \(\ell\)
one has
\[
 \operatorname{ord}_\ell(d)\equiv
       -\operatorname{ord}_\ell(Y_0)\pmod P,\qquad
 \operatorname{ord}_\ell(d)\equiv
       -\operatorname{ord}_\ell(Z_0)\pmod Q .
\]
The moduli are coprime. At all but finitely many primes
the prescribed class modulo \(PQ\) is zero. Choosing
the least nonnegative solution at every prime gives
a fixed \(d_0\ge1\) such that every admissible \(d\)
has shape \(d_0 t^{PQ}\) for a positive integer \(t\).
The divisibility \(d^g\mid j\) now gives
\(t\ll X^{s/(gPQ)}=X^{s/\operatorname{lcm}(p,q)}\).
The lifting maps are injective on positive integer roots,
so this is the desired number of original root pairs.
The same combination of two energies gives the
seven-arm affine count, and taking the largest
admissible \(p,q\) is strongest for fixed groups.

#### An exact different-power example

Consider the unit-budget profile
\[
 a_1=11/100,\quad a_2=3/50,\quad a_4=7/50;\qquad
 b_2=17/100,\quad b_3=1/50,\quad b_5=3/50,\quad b_6=1/20;\qquad
 c_1=7/200,\quad c_2=27/100,\quad c_5=17/200.
\]
The masses of the letters are \(31/100,3/10,39/100\).
The last letters both have weighted height one, and
the first has height \(79/100\).
Lock \(b_5,b_6\), giving \(f=11/100,w=3/5,s=2/5\).

For \(a\), use \(I=\{2,4\}\) with \(p=2\) and
\(J=\{4\}\) with \(q=4\). The lifted root heights are
\(34/100,28/100\), so the rank-two gap is \(2/100\).
The collinear gap is \(1/10\), and (B) gives
\(E_a=19/50\). The best equal-power choice, allowing
asymmetric groups too, has \(E_a=21/50\) and uses
\(\{2,4\}\) in both copies with power 2. This follows
by checking the three permitted groups \(\{2\}\),
\(\{4\}\), and \(\{2,4\}\).
On \(c\) choose \(\{2\}\) in both copies, giving
\(E_c=11/25\). This is also the best equal-power choice
on that side.

The resulting cap is \(13/25\), compared to \(27/50\)
with the best equal-power choices for the same lock.
This is a comparison of the two repeat bounds at a
fixed modulus.

### Slicing the rank-two lattice for unequal lifted powers

In the genuinely different-power repeat bound, the two copies
release \(p\)-th and \(q\)-th power factors of root exponents
\(r_1,r_2\). Write \(g=\gcd(p,q)\), \(P=p/g\), \(Q=q/g\)
and lift to \(Y=u^P,Z=v^Q\). The raw lifted exponents are
\(\rho_1=P r_1,\rho_2=Q r_2\).
Suppose the lock supplies a modulus of numerical exponent
at least \(w>\max(\rho_1,\rho_2)\).
Write \(s\) for the upper exponent of the nonzero repeat
quotient and \(t=\rho_1+\rho_2-w\).
The repeat bound can use the gap
\[
 \min\left\{\max\{t,s/\operatorname{lcm}(p,q)\},
  \max\left\{
        t/2+\min(r_1/Q,r_2/P),
        s/\operatorname{lcm}(p,q)\right\}\right\}.
 \tag{S}
\]
As usual add \(2L_T-m_I-m_J\) to obtain the repeat energy
exponent. Formula (S) is a choice between valid maxima
of affine terms when the two group powers are fixed.

For fixed complementary coefficients, the unit ratio congruence
for \(Y/Z\) is a fixed \(g\)-th power modulo the lock.
Thus all root pairs fall in \(O_\varepsilon(X^\varepsilon)\)
lattices of index \(M\), each of the form
\(Y\equiv\lambda Z\bmod M\) with \(\lambda\) a unit.
Normalize the lifted box of bounds \(U,V\) to a fixed
square. With two independent vectors in the square,
let \(\ell\) be the shortest vector length and \(H\) the
perpendicular component of a reduced second basis vector.
Both are bounded above, \(H\gg\ell\), and
\(\ell H=M/(UV)\). The coefficient of the second basis
vector is bounded by \(O(1/H)\ll (UV/M)^{1/2}\),
so this coordinate has \(O(X^{t/2})\) possible values.
When \(t<0\) strictly, this two-dimensional case cannot
occur for large \(X\); at a vanishing exponent slack
the same estimates allow \(X^{o(1)}\) errors.

The shortest vector is nonaxial. In original coordinates
a nonzero horizontal or vertical lattice vector has its
nonzero coordinate at least \(M\), since \(\lambda\)
is a unit. Its length after normalization would be
at least \(M/U\) or \(M/V\), respectively, which tends
to infinity by the strict guard. Thus fixing a nonzero
value of the second integer coefficient cuts the lifted
points by an affine line
\(aY+bZ=c\) where \(a,b,c\ne0\). In the original
roots this is \(a u^P+b v^Q=c\).
It is absolutely irreducible, by Eisenstein at a simple
root of the binomial in the other variable. The lattice
basis in original coordinates and the line constant have
polynomial height in \(X\) because both normalized basis
vectors have bounded length. The plane rectangle bound using the two pure monomials gives
\[
 O_\varepsilon(X^{\min(r_1/Q,r_2/P)+\varepsilon})
\]
points on this line pullback. The line through the origin
has at most the collinear repeat bound, since its primitive
ray forces the scale to be a fixed coefficient times a
\(PQ\)-th power and the scale to the \(g\)-th power
divides the nonzero repeat quotient. The same collinear
proof covers the whole box if it contains no independent
vectors. This proves the second alternative in (S).
The first alternative is the earlier unsliced lattice
bound.

The proof also applies to primitive cross-letter lattice
pairs, where the zero line has at most one pair, so one
may replace \(s/\operatorname{lcm}(p,q)\) by zero.
For the short-addend real strip on such a pair, use
the real linear change with coordinates
\[
 (Y-\kappa Z)/X^{\rho_1+h-1},\qquad Z/X^{\rho_2},
\]
where \(h\) is the short height and \(\kappa\asymp
X^{\rho_1-\rho_2}\) comes from equality of the tall
powers. The lock is the short cofactor, of weight
\(h-k r\) if the released short term is a \(k\)-th
power with root exponent \(r\). The normalized covolume
gives \(t=\rho_1+\rho_2+kr-1\) and an axis vector would
have length at least a constant times
\(X^{1-kr-\rho_1}\) or \(X^{1-kr-\rho_2}\). Thus
the projected nonaxial guard is \(kr+\max(\rho_1,\rho_2)<1\).
Both formulas cancel the unknown short tail weight.

#### Exact unit-budget comparison

Take \(a_1=53/750,a_5=19/250,a_6=19/300\), so
\(L_a=21/100\) and \(a\) is short. Take
\(b_1=94/1000,b_2=103/1000,b_6=21/1000,b_7=82/1000\),
so \(L_b=3/10\) and \(b\) has height one. Take
\(c_1=1127/8000,c_2=591/2000,c_5=429/8000\), so
\(L_c=49/100\) and \(c\) also has height one. Lock \(b_6\)
and the whole tail of \(b\), with
\(f=103/1000,w=7/10,s=3/10\). The total tail weight is determined by the retained weights
on the tall letter.

Release the degree-5 block of the first copy of \(a\)
and the degree-6 block of the second copy at their full
powers. The lifted heights are both \(19/50\) and the
axis guard is strict. Then \(t=3/50\), and (S) takes
the sliced rank-two arm \(16/375\) over the collinear
arm \(1/100\). This gives \(E_a=97/300\). The best
unsliced energy on \(a\), allowing different groups
and powers including first powers, is \(41/125\).
Both copies using degree five attain this, as does
using \(\{1,5\}\) with power one in the first copy
and degree five in the other. On \(c\), both copies
use the degree-2 block, giving the optimal
\(E_c=539/1000\) with or without slicing. Use
\(A=c,B=a\) in the separated-diagonal bound to
obtain
\[
 D\le641/1200<107/200,
\]
compared with \(1073/2000>107/200\) for the best
unsliced cap with this modulus.

[Heath-Brown's affine rectangle estimate, Equation (3) in the joint paper](https://arxiv.org/pdf/math/0506497)

### Exact valuation support in a mixed residue ray

Use the anatomic coordinates \(T=T_0\prod_i T_i^i\), with
the bases \(T_i\) pairwise coprime and with the early indices
\(1,\ldots,6\) exact. In the two copies of an off-diagonal
repeat on \(T\), release valuation groups \(I,J\) at fixed
powers \(p,q\), respectively. An entire aggregate tail, when
released, is assigned power one and represents the indices
\(7,\ldots,d\). Write \(g=\gcd(p,q)\), \(P=p/g,Q=q/g\),
and lift the two roots to \(Y=u^P,Z=v^Q\).
An index \(i\) in the first group contributes a prime exponent
at least \(i/g\) on the lifted value, with equality if its base
is squarefree, and similarly for \(J\).
Let \(m_I,m_J\) be the two released masses. Let \(w\) be a
lower numerical weight for the modulo lock \(M\), and let
\(s\) bound the logarithmic magnitude of the nonzero quotient
\(j\) in the complete repeat difference \(jM\).

Let \(e\) be the least possible shared valuation index:
\[
 e=\min(I\cap J),\qquad s/e:=0\quad\hbox{if }I\cap J=\varnothing .
\]
When the only overlap is an aggregate tail, use \(e=7\).
The old collinear arm \(s/\operatorname{lcm}(p,q)\)
can be replaced by \(s/e\), without altering the rank-two arms.
In particular this gives
\[
 E_T\le 2L_T-m_I-m_J+
     \max\{\rho_1+\rho_2-w,\ s/e\},
\]
where \(\rho_1=P r_I,\rho_2=Q r_J\) are the lifted root
exponents. A real-strip or sliced rank-two bound can be
substituted on the same terms. Since any shared index is a
multiple of \(\operatorname{lcm}(p,q)\), this never worsens
the old collinear arm.

#### Proof and exceptional coefficient

Take the constructive BBLT Lemma 2.2 decomposition with
cutoff parameter \(k=\lceil\epsilon_0^{-1}\rceil>6\).
Its only nonsquarefree base can be \(T_k\), which absorbs
integer powers from valuations above its large cutoff.
Thus the retained bases of indices \(1,\ldots,6\) are
squarefree and genuinely of their stated indices.
All hidden coordinate exponents are at least seven.
The exceptional factor \(T_0\) is fixed with the unreleased
coordinates and costs \(X^{O(\epsilon_0)}\), which is
absorbed as usual. Its prime support can also be included
in the fixed support below, including the exceptional
factor of each copy if they are not identical.
The statement only uses exact early valuations outside
this fixed exceptional support. [BBLT v2, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

Freeze the two complementary products. Unit ratio classes
of \(Y/Z\) modulo \(M\) give \(X^\epsilon\) lattices of
index \(M\), as in the mixed-power bound.
In a collinear lattice write \((Y,Z)=d(Y_0,Z_0)\)
with \((Y_0,Z_0)\) its positive primitive integer ray.
The fixed coefficient equation gives \(d^g\mid j\)
because \(\gcd(d,M)=1\) and the repeat is nonzero.
Let \(R\) be the product of primes dividing \(Y_0Z_0\)
or either exceptional coefficient. It is at most \(X^{O(1)}\).
Write \(d=d_R d_*\) according to this fixed support.

For a prime \(\ell\nmid R\) dividing \(d_*\), both
lifted powers have valuation \(v_\ell(d_*)>0\). It
follows that their released numeric products have the
same complete exponent \(g\,v_\ell(d_*)\). If one of
the copies puts this prime in a retained early coordinate,
the other copy must put it in the same retained coordinate:
a hidden coordinate would give exponent at least seven.
If both put it in hidden coordinates, its exponent is
at least seven. Thus the exponent is at least \(e\).
When the two groups have no common retained index and
do not both contain the tail, \(d_*=1\).

When \(e<\infty\), we have
\[
 \operatorname{rad}(d_*)^e\le d_*^g\le |j|\ll X^s .
\]
There are \(O(X^{s/e})\) radical choices. Both the
fixed-support multiplicity \(d_R\) and the multiplicity
of \(d_*\) at a fixed radical cost \(X^\epsilon\):
for any \(R'\le X^C\) and \(\tau>0\),
\[
 \#\{m\le X^C:\operatorname{rad}(m)\mid \operatorname{rad}(R')\}
 \le X^{C\tau}\prod_{\ell\mid R'}(1-\ell^{-\tau})^{-1}
 \ll_\tau X^{2C\tau}.
\]
For the final inequality use
\((1-\ell^{-\tau})^{-1}\le \ell^\tau\) for all sufficiently
large primes, and absorb the finitely many remaining
factors in the constant. Taking \(\tau\) sufficiently
small gives the claim, uniformly in the ray and fixed
exceptional coefficients. Summing over roots of the
unit congruence and the two frozen complementary
products proves the refined collinear arm.

#### Exact fixed-lock comparison

Take this unit-budget profile:
\[
 a_1=1/25,\quad a_2=7/50,\quad a_4=4/25;\qquad
 b_1=1/250,\quad b_4=b_5=1/10,\quad b_6=2/125;\qquad
 c_1=21/800,\quad c_2=73/200,\quad c_5=39/800.
\]
The masses are \(17/50,11/50,11/25\) and both \(b,c\)
have height one. The first letter has height \(24/25\).
Lock \(b_5\), yielding \(f=1/10,w=s=1/2\).
For \(a\), take disjoint groups \(\{2\},\{4\}\), at powers
two and four. The lifted exponents are \(7/50,8/25\),
whose sum \(23/50\) is strictly less than \(w\).
The rank-two arm is negative and the collinear arm is
zero, giving \(E_a=19/50\). For \(c\), use the degree-2
block in both copies, giving \(E_c=2/5\).
With \(A=c,B=a\) the count exponent is at most \(49/100\).

Enumerating all nonempty exact groups and all allowed
powers on the two copies, including first powers and
the guarded sliced arm, gives \(E_a\ge97/200\)
and \(E_c\ge2/5\) with the old collinear arm. Optimizing
the separated-diagonal orientation then gives
\(217/400>107/200\) with the old bound.
The comparison fixes the modulus.

### Sliced primitive ratio locks with axial lines

Release groups on the free letters \(A,B\) in powers \(p,q\),
and let \(r_1,r_2\) be their original root exponents. Choose
any common divisor \(g\) of \(p,q\), put \(P=p/g\),
\(Q=q/g\), and lift the roots to \(Y=x^P,Z=y^Q\) with
exponents \(\rho_1=P r_1,\rho_2=Q r_2\).
Freeze the two complementary factors and the chosen lock
on the third letter. The log cost of these choices is
\(F=L_A+L_B-m_A-m_B+m_{\rm lock}\).

For a modulus \(M\) of lower weight \(w\), let
\(u=w\). If the locked letter is the short addend and its
complement has upper numerical weight \(s\), we may instead
take \(u=1-s\). In either case put
\(t=\rho_1+\rho_2-u\). Then the ordinary and sliced bounds
are respectively
\[
 D\le F+(t)_+,\qquad
 D\le F+
 \left[\frac t2+\max\!\left\{
   \min\!\left(\frac{r_1}{Q},\frac{r_2}{P}\right),
   r_1-u,\ r_2-u
 \right\}\right]_+ .
 \tag{RS}
\]
The maximum can use either of the two terms in the inner
minimum as an affine upper bound. Both estimates allow
arbitrarily small exponent losses.

As in the mixed-power residue count, reduction modulo \(M\)
gives \(O_\varepsilon(X^\varepsilon)\) unit ratio lattices
of index \(M\).
For the ordinary estimate, normalize the box on \(Y,Z\)
to bounded side lengths. If the lock is short, the actual
two tall values are comparable to \(X\) and differ by at
most \(O(MX^s)\). For fixed coefficients their normalized
lifted pair lies in a parallelogram specified by
\[
  |Y-\kappa Z|\ll X^{\rho_1}MX^{s-1},\qquad
  Z\ll X^{\rho_2}, \qquad \kappa\asymp X^{\rho_1-\rho_2}.
\]
This follows by taking the ratio of the two positive
\(g\)-th powers and using the mean value theorem on a
fixed compact interval. Upon scaling this parallelogram
to bounded size the lattice has determinant
\(\gg X^{-t}\), since its index \(M\) cancels the \(M\)
in the strip width. The ordinary normalized box also has
determinant \(\gg X^{-t}\).

If the normalized region contains two independent vectors,
use a Gauss-reduced basis with shortest vector of length
\(\ell\) and perpendicular spacing \(H\) between its lines.
Both are bounded above, \(H\gg\ell\), and
\(\ell H\) is the lattice determinant. The area estimate
counts \(O((\ell H)^{-1})\) points. Alternatively, there
are \(O(H^{-1})\ll X^{t/2}\) parallel lines. If \(t<0\)
strictly, two independent vectors are excluded for
large \(X\). In the collinear case the lifted pair is
coprime because the two complete letters are coprime.
It therefore lies at the unique positive primitive
point on its ray.

For the sliced estimate, the line through the origin
again contributes at most one coprime pair. On a
nonaxial line not through the origin the equation for
the original roots is \(\alpha x^P+\beta y^Q=\gamma\)
with all three coefficients nonzero. Its irreducibility
follows by Eisenstein at a simple zero of
\(\beta y^Q-\gamma\), including the case \(P=1\).
The affine rectangle count bounds the number on
one line by \(O_\varepsilon(X^{\min(r_1/Q,r_2/P)+\varepsilon})\).
[Browning–Heath-Brown, equation (3)](https://arxiv.org/pdf/math/0506497)

If the line is axial, one lifted root is fixed and the
other is in a fixed unit residue class modulo \(M\).
For instance, at fixed \(Z\), there are \(O_\varepsilon(X^\varepsilon)\)
possible \(x\) residue classes modulo \(M\). This counts
\(O_\varepsilon(X^{(r_1-w)_++\varepsilon})\) choices of \(x\)
in the ordinary box. In the short case the difference of
the complete values bounds the interval of \(x\) by
\(O(X^{r_1} M X^{s-1})\), by the mean value theorem
on the original tall powers. Thus in this case there are
\(O_\varepsilon(X^{(r_1+s-1)_++\varepsilon})\) choices.
The same arguments apply with \(x,y\) interchanged.
This accounts for both axis terms in (RS), with no
separate restriction on which direction is shortest.

For an exact comparison, take the nonzero unit-budget
masses
\[
 a_1=1/25,\quad a_6=3/20,\qquad
 b_1=11/200,\quad b_2=13/50,\quad b_5=17/200,\qquad
 c_1=1/50,\quad c_2=29/100,\quad c_4=1/10 .
\]
The short height is \(47/50\) and both tall heights are one.
Lock \(c_4\) and release \(a_6,b_2\) with common divisor \(2\).
Here \(w=2/5,\rho_1=9/20,\rho_2=13/50\), \(F=7/25\),
and \(t=31/100\). The ordinary ratio bound with these
groups is \(59/100\). The sliced bound has curve term
\(13/150\) and both axis terms are negative, giving
\(313/600<107/200\). In this example \(w<\rho_1\).

For the short-addend form, take
\(a_1=a_3=b_4=c_4=1/10\) and \(b_2=c_2=3/10\).
The \(a\)-height is \(2/5\) and the other two heights are one.
Lock \(a_3\) and release \(b_2,c_2\) with \(g=2\).
The ordinary margin is \(3/10\), whereas the short margin
is \(-3/10\), so (RS) gives \(3/10\) instead of
\(3/5\) for this fixed choice of groups and lock.

## Signed monomial congruences over full moduli

Fix \(m\ge1\) and \(C>0\). For all sufficiently large \(k\), depending only
on \(m,C\), let \(n\ge2\) be \(k\)-full, so \(p^e\Vert n\) implies
\(e\ge k\). Fix nonzero signed integer exponents \(\mathbf a=(a_i)_{i=1}^m\), and let \(\mathcal S\) be a set of unit tuples of positive integers
\(x_i\le n^{C/k}\) such that
\[
\phi(\mathbf x)=\prod_i x_i^{a_i}\equiv K\pmod n.
\]
Suppose the exact rational monomial fibers in \(\mathcal S\) have size at most
\(F_0\). Then, uniformly in the unit target,
\[
\#\mathcal S\ll_{m,C,k,\mathbf a,\varepsilon}
F_0 n^{m(m-1)/(2k)+\varepsilon}.
\tag{FC1}
\]
The lower bound on \(k\) is independent of \(\mathbf a\).

The count uses an induction on algebraic dimension. The logarithmic Jacobian
controls the amount of lost precision across the prime powers of \(n\). The
coefficients of this logarithmic Jacobian admit a height bound uniform in the
target.

### The logarithmic Jacobian and local lifts

Start with \(Z=\mathbb A^m\) and work at a stage on an algebraic set \(Z\) of
dimension \(j\). A zero-dimensional set has boundedly many points. At positive
dimension, the integer defining equations have bounded degrees and heights
\(n^{O(1/k)}\). Choose a projection \(\mathbf t\) of \(j\) bounded integer
linear forms dominant on every top-dimensional component, completing it to an
invertible bounded integer linear change. The choice is by nonvanishing of
finitely many tangent minors and the determinant, whose total degree is
bounded. Elimination gives a nonzero equation in \((\mathbf t,z_i)\) for each
of the \(m-j\) trailing coordinates, since each such projection has dimension
at most \(j\). Remove the factors in \(\mathbf t\) alone and take the
squarefree part \(G_i\) in \(z_i\). The removed factors and the derivatives
\(\partial G_i/\partial z_i\) vanish only on proper subsets of each top
component, by dominance of \(\mathbf t\) and squarefreeness over
\(\mathbb Q(\mathbf t)\). Count these lower-dimensional sets by induction.
Pull the \(G_i\) back to the original coordinates by the integer change.
Hensel lifting then uses only integer coefficients, including at primes
dividing the determinant of the linear change.

On the remaining open set the differentials of the \(G_i\) are independent.
Points on components with constant \(\phi\) cost at most \(O(F_0)\). On the
other components take a bounded-integer linear combination \(J\) of maximal
minors of the rows
\[
\nabla G_1,\ldots,\nabla G_{m-j},
\quad
\left(a_i\prod_{h\ne i}X_h\right)_{i=1}^m
\]
that is generically nonzero on each remaining component. Indeed, if the rank
is deficient along such a component, the logarithmic differential annihilates
its tangent space, so \(\phi\) is constant there in characteristic zero.
Clearing denominators in that differential leaves a row of degree \(m-1\).
The \(J=0\) points on the nonconstant components lie in the intersection of
\(J=0\) with the closure of \(Z\setminus\{J=0\}\), which has dimension at
most \(j-1\). Elimination and saturation preserve bounded degree and height
\(n^{O(1/k)}\). These proper exceptional sets are handled by induction.

The coefficient-height bound for elimination follows by stratifying the input
coefficient space into finitely many locally closed pieces on which a Gröbner
basis specializes from fixed rational functions in the inputs. Compute
generically over each parameter component's function field, exclude the
finitely many denominators in the basis and its membership and reduction
identities, and continue on the closed complement. Clearing denominators gives
an exponent \(O(1/k)\) for the output integer coefficients. Choosing a
primitive factor or a squarefree part preserves that exponent by the
fixed-degree factor-height bound, for example by Kronecker substitution
followed by its univariate form. Saturation by a polynomial \(J\) is obtained
by adding \(wJ-1\) and eliminating \(w\). These operations bound the degrees
and heights in terms of the algebraic input bounds. The logarithmic row uses
the fixed integers \(a_i\), whose sizes can be absorbed into a factor
\(n^{1/k}\) for sufficiently large \(n\).

At a remaining integer point \(\mathbf x\) write \(|J(\mathbf x)|\le n^{A/k}\), with \(A\) a constant depending only on the input bounds, \(m,C\), after
absorbing fixed constants for large \(n\). For \(p^e\Vert n\), record \(p\)
as bad if \(v_p(J(\mathbf x))>e/3\). The discarded prime powers have product
at most \(n^{3A/k}\), because
\[
\sum_{p\mid n}v_p(J(\mathbf x))\log p\le(A/k)\log n.
\]
At any good prime at least one of the chosen minors has valuation at most
\(e/3\). Choose an integer representative for \(K\) and set
\(F=\prod X_i^{a_i^+}-K\prod X_i^{a_i^-}\), where \(a_i^+=\max(a_i,0)\) and
\(a_i^-=\max(-a_i,0)\). At the original unit point, \(\nabla F\bmod p^e\) is
a unit multiple of the cleared logarithmic row. The valuation \(\delta\) of
the corresponding minor of \((\nabla G_i,\nabla F)\) is therefore at most
\(e/3\). Keep the unselected original coordinates fixed and lift the system
\(G_i=F=0\) by multivariate Hensel to a unit point of distance at most
\(p^{-(e-\delta)}\). Here \(e>2\delta\). The lifted minor stays nonzero. The
reference set given by \(G_i=F=0\) and the full-rank condition therefore has
dimension \(j-1\) at each of its points.

### The determinant count and amplification

For \(j>1\), the exact reference sets form a semialgebraic family of dimension
\(j-1\) when nonempty. For sufficiently large primes, uniform \(T_r\)
parametrization gives boundedly many charts for every fixed order \(r\),
uniformly in the defining coefficients. Over each of the finitely many smaller
fields \(\mathbb Q_p\), the fixed-field parametrization theorem indexes its
charts by a fixed finite power of a finite residue-ring sort. Its chart count
is therefore also uniform in the coefficient parameters. The small-prime chart
constants contribute at most a fixed product, independently of their powers in
\(n\).
[Cluckers–Forey–Loeser, Theorem 3.1.4](https://arxiv.org/pdf/1902.06589v1#page=15)
·
[Cluckers–Comte–Loeser, Theorem 3.2.3](https://arxiv.org/pdf/1404.1952v3#page=24)

A \(T_r\)-chart restricted to a parameter class modulo \(p^b\) agrees modulo
\(p^{rb}\) with a polynomial of degree less than \(r\) and integral
coefficients in the parameter differences. This is the approximation used
below when the original points are replaced by their corrected points.

At each good prime partition parameter charts modulo \(p^{b_p}\), where
\(b_p=\lfloor e/k\rfloor\ge1\). Record the chart and parameter class of the
chosen correction, or that the prime was bad. The total number of records is
\[
O_\varepsilon\left(n^{(j-1)/k+\varepsilon}\right).
\]
Indeed each large prime contributes at most \(1+C_r p^{b_p(j-1)}\), the
product of the \(p^{b_p}\) is at most \(n^{1/k}\), and
\(C_r^{\omega(n)}\ll_\varepsilon n^\varepsilon\). Small primes have the
separately bounded constants just described. In a given record the set of
ignored primes is the same for all points.

Use the \(D=\binom{M+j}{j}\) monomials of total degree at most \(M\) in the
\(j\) projected coordinates \(\mathbf t\). Write \(w_1,\ldots,w_D\) for the
least total degrees of \(D\) distinct monomials in \(j-1\) Taylor variables
and put \(S_D=\sum w_i\). Since \(S_D\asymp_j D^{1+1/(j-1)}\), choose \(M\)
so that \(S_D>4CMD+4\), and take \(r>w_D\). Increasing \(k\) in terms of
these constants ensures \(2e/3\ge r b_p\) and \(1-3A/k\ge1/2\). An evaluation
determinant from one record is divisible at good \(p\) by \(p^{b_p S_D}\).
Taylor expansion and multilinearity give this divisibility. Repeated Taylor
monomial columns give zero, while every error column, including the lifting
error, is divisible by \(p^{rb_p}\). The projected coordinate forms have
integer coefficients, so substitution and taking monomials preserve the
congruent Taylor approximations. Since \(b_p\ge e/(2k)\), the simultaneous
divisibility is at least
\[
n^{(1-3A/k)S_D/(2k)}\ge n^{S_D/(4k)}.
\]
The determinant is an integer of absolute value \(O(n^{CMD/k})\), so
vanishes. Its nullspace gives a bounded-degree polynomial in \(\mathbf t\) of
height \(n^{O(1/k)}\) by maximal minors. It cuts dimension on every top
component by dominance of the projection. This produces
\(O_\varepsilon(n^{(j-1)/k+\varepsilon})\) child sets of dimension at most
\(j-1\).

At \(j=1\), the good-prime corrected points are isolated solutions of
\(G_i=F=0\). Their number is bounded uniformly in the coefficients by Bézout,
with the bound allowed to depend on the fixed exponents \(\mathbf a\). Record
the isolated solution or a bad prime at each \(p\mid n\). There are
\(O_\varepsilon(n^\varepsilon)\) possibilities. Two original points with the
same record agree modulo \(p^{\lfloor 2e/3\rfloor}\) at every good prime. The
combined modulus is larger than \(2n^{C/k}\) for large \(k\), because the
total good logarithmic mass is at least \((1-3A/k)\log n\) and
\(\lfloor2e/3\rfloor\ge e/2\). Hence each record fixes at most one original
tuple.

Summing the dimension costs gives \(\sum_{j=1}^m(j-1)=m(m-1)/2\), proving
(FC1). The degree and height bounds, and the recursively required minimum
\(k\), depend only on \(m,C\), because all logarithmic-minor and auxiliary
degrees are independent of \(\mathbf a\). The large-prime cutoff and
fixed-prime chart constants may depend on \(\mathbf a\) and the chart order
and only enter the implied constant.

For integer tuples with pairwise coprime coordinates, group solutions by a
prime coloring and amplify by taking coordinatewise products. Assign one color
to each coordinate. Some coloring retains \(n^{-o(1)}\) of the tuples with
each prime in its coordinate's color. Form \(\ell\)-fold products from those
tuples. Input multiplicity is divisor-many, and the signed exact monomial on
the colored products has fibers of size one, because the valuation of a prime
is assigned to its unique coordinate. Apply (FC1) to the products at height
constant \(\ell C\) and target \(K^\ell\), then take \(\ell\)-th roots. The
bound becomes
\[
O_{m,C,k,\mathbf a,\ell,\varepsilon}
\bigl(n^{m(m-1)/(2\ell k)+\varepsilon}\bigr)
\tag{FC2}
\]
for sufficiently large \(k\) depending on \(m,\ell C\). Thus any bounded
collection of pairwise coprime coordinates of height \(X^{O(L)}\) can be
released at cost \(o(L)\) when the modulus has positive limiting weight and
fullness threshold of order \(1/L\). A grouped low cofactor may be treated as
a coordinate of exponent one when its height is \(X^{O(L)}\). In the anatomic
decomposition we fix the small exceptional coefficient at the truncation loss.
The remaining coordinates within each letter may be taken pairwise coprime.
Coordinates on different letters are coprime by primitivity. All exponents
range over a fixed finite set at each truncation before \(X\to\infty\).
[Bernert–Browning–Lichtman–Teräväinen, Proposition 2.1 and Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

## Signed monomials in real short strips

Fix \(m\ge1\) and \(\eta>0\). There exists \(K_0=K_0(m,\eta)\) with the
following property. Fix any signed integer vector
\(\mathbf a=(a_1,\ldots,a_m)\), set \(\phi(\mathbf x)=\prod_i x_i^{a_i}\),
and let \(C>0\) vary arbitrarily. For \(K\ge K_0\), consider positive integer
tuples \(x_i\le H\) such that
\[
|\phi(\mathbf x)/C-1|\le H^{-K}.
\tag{SR2}
\]
If every exact rational value of \(\phi\) has at most \(F_0\) preimages in the
set under consideration, the count is
\[
O_{m,\eta,\mathbf a,K}(F_0 H^\eta),
\tag{SR3}
\]
uniformly in \(C\). In particular \(K_0\) is independent of the signed
exponent vector.

On each dyadic box normalize \(z_i=x_i/B_i\in[1,2]\), with integer
\(B_i\le H\). Put
\[
\Phi(\mathbf z)=\frac{\prod_i B_i^{a_i}}{C}\prod_i z_i^{a_i}.
\]
The scalar factor is bounded above and below in terms of the fixed exponents
if the box has at least one solution. Work on the enlarged positive box
\([1/2,3]^m\). The exact reference equation is \(\Phi=1\), a polynomial
equation after clearing powers of the coordinates. Its degree may depend on
\(\mathbf a\), but its logarithmic row has the form
\[
d\log\Phi=\left(a_i/z_i\right)_{i=1}^m,
\]
which clears to a polynomial row of degree \(m-1\).

Apply the algebraic dimension reduction of the full-modulus proof to the
bounded-degree sets \(Z\) containing the normalized rational points, starting
from \(Z=\mathbb A^m\). Their defining coefficients have height \(H^{O(1)}\),
with the exponent independent of \(\mathbf a,C\). At a set of dimension \(j\),
choose \(j\) bounded integer linear projection forms \(\mathbf t\) dominant on
each top-dimensional component, as in (FC1). Eliminate to get polynomial
equations \(G_i(\mathbf t,y_i)\) for the \(m-j\) trailing coordinates, with
bounded degree and height \(H^{O(1)}\). The uniform Gröbner-basis
specialization above gives the coefficient-height bound. Remove the
\(\mathbf t\)-only factors and take the squarefree part in \(y_i\). The zero
sets of the removed factors and the trailing derivatives \((G_i)_{y_i}\) meet
\(Z\) in sets of dimension at most \(j-1\) by dominance and squarefreeness, so
they are counted by induction. On the remaining set, the \(G_i\), pulled back
to the normalized variables, have independent gradients.

An irreducible geometric component on which \(\Phi\) is constant contributes
at most \(F_0\) points. On any other component meeting the good open set, the
cleared logarithmic row together with the \(\nabla G_i\) has a generically
nonzero maximal minor. Otherwise the logarithmic differential annihilates that
component's tangent space and the monomial is constant, in characteristic
zero. Choose a bounded-integer linear combination \(J\) of the minors that is
nonzero generically on each such component. The remaining bad \(J=0\) points
are in the intersection of \(J=0\) with the saturation of \(Z\) by \(J\), a
bounded-degree set of dimension at most \(j-1\) and with coefficient height
\(H^{O(1)}\). Use induction there as well.

At every remaining original point, some minor of \((\nabla G_i,d\log\Phi)\) is
nonzero, with inverse absolute value \(O(H^A)\) for a fixed \(A\) independent
of \(\mathbf a\). Indeed each \(z_i\) is rational with denominator at most
\(H\). After clearing the factors from the logarithmic row, a nonzero
integer-polynomial minor of bounded degree has absolute value at least
\(H^{-O(1)}\) at such a point. Dividing the \(G_i\) by their coefficient
maxima to keep their derivatives bounded loses only another \(H^{O(1)}\). The
linear change used in the projection has bounded coefficients. The exponents
\(a_i\) contribute only fixed constants to derivative bounds, not to their
powers of \(H\).

Fix the remaining \(j-1\) original coordinates and solve the real system
\[
G_i=0,\qquad \log\Phi=0
\]
in the \(m-j+1\) coordinates specified by the minor. At the original point its
residual is \(O(H^{-K})\) and the inverse Jacobian has norm \(O(H^A)\). Its
second derivatives are bounded on the enlarged box, with constants allowed to
depend on the fixed exponents. For \(K>2A+1\), quantitative inversion gives a
solution within \(O(H^{-K+A})\) of the original point, and the relevant
Jacobian remains nonzero there. This follows directly from the contraction
mapping theorem applied in a ball of radius a sufficiently large fixed
multiple of \(H^{-K+A}\), using the original inverse Jacobian. The exact
solution belongs to the semialgebraic set defined by the \(G_i=0\), the
cleared equation \(\Phi=1\), and the full rank condition. That reference set
has real dimension \(j-1\) on the enlarged positive box. Use this exact
reference set for the determinant comparison.

When \(j>1\), the real Yomdin–Gromov algebraic lemma gives boundedly many
charts of dimension \(j-1\) with bounded derivatives through any fixed order.
The bounds are uniform in the coefficients, with the semialgebraic complexity
fixed by the equation degrees.
[Binyamini–Novikov, Theorem 1](https://arxiv.org/pdf/2009.02082v1#page=1)

Divide the parameter domains into cells of side \(H^{-\tau}\), with fixed
\(\tau>0\). Use the \(D=\binom{M+j}{j}\) monomials of total degree at most
\(M\) in \(\mathbf t\). If \(b_1,\ldots,b_D\) are the first Taylor orders in
\(j-1\) variables, the determinant in any one cell is
\[
O\bigl(H^{-\tau\sum_{\nu=1}^D b_\nu}\bigr)
\]
as long as \(K-A\ge\tau b_D\). Taylor-expand about a parameter point in that
cell. Repeated Taylor monomial vectors give zero contributions and the error
in each column has at least the weight of the largest order used. The
correction error of size \(O(H^{-K+A})\) is allowed by the same inequality.

Each projected coordinate at an original point has denominator dividing
\(\prod_i B_i\le H^m\). Clearing the determinant denominators costs at most
\(H^{mMD}\). Choose \(M\) so that \(\tau \sum b_\nu>mMD\), possible because
the number of projected coordinates exceeds the parameter dimension. The
integer determinant is then zero. A maximal minor of the original evaluations,
with denominators cleared, gives a nonzero polynomial in \(\mathbf t\) of
bounded degree and height \(H^{O(1)}\). Its vanishing cuts all top-dimensional
components of \(Z\) by dominance of the projection. There are
\(O(H^{\tau(j-1)})\) such child sets, with constants allowed to depend on the
fixed exponent vector. This induction, including the finite additional
lower-dimensional exceptional sets, incurs at most \(H^{\tau m(m-1)/2}\) in
total.

When \(j=1\), the corrected points are isolated algebraic solutions of
\(G_i=0,\Phi=1\) by the full-rank Jacobian. Bézout bounds their number in
terms of the equation degrees, including the fixed exponents. For \(K>A+2\)
their correction neighborhoods contain at most one original point each, as
distinct normalized rational points from the same dyadic box are separated by
at least \(1/H\). Dimension zero needs no correction. Choose \(\tau\)
sufficiently small in terms of \(m,\eta\) and absorb the dyadic count into the
remaining exponent. In the induction the determinant degrees, height exponents
and the real correction precisions needed depend only on \(m,\eta\). The
potentially high degree of the reference equation affects only the chart
constants and the terminal Bézout bound. This proves (SR3).

In particular, when the short term has limiting weight zero, the two tall
terms have ratio \(1+O(X^{-1+o(1)})\). Suppose both tall terms have
bounded-coordinate representations with all variable coordinates of height
\(X^{O(L)}\), with any fixed coefficient factors allowed. Divide the two terms
to get (SR2) with signed exponents. The target \(C\) absorbs the fixed
factors, and the precision tends to infinity as \(\lambda\to0\), taking
\(H=X^{O(L)}\). In each exact fiber the ratio of the two complete coprime
terms determines both terms. Their anatomic representations are divisor-many
at fixed truncation. Thus the two tall terms can be counted together at
normalized cost zero after any necessary factors have been fixed.

## A limiting small-budget coefficient of \(1/3\)

Let \(N_\lambda(X)\) count ordered positive primitive triples \(a+b=c\le X\)
with \(\operatorname{rad}(abc)<c^\lambda\). For every \(\eta>0\) and all
sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\eta,\varepsilon}
X^{(1/3+\eta)\lambda+\varepsilon}.
\]
In the limiting box notation \(L\) is the total root mass, \(l_T\) is the
normalized mass on \(T\), and \(d=D/L\) in the limit. The numerical weight of
a factor is its logarithmic exponent to base \(X\). Suppose \(d>1/3\), put
\(u=1-2d\), and call a letter supported when \(l_T>u\). The fixed-degree
cutoff gives a supported letter boundedly many heavy roots and a grouped low
cofactor, all of height \(X^{O(L)}\), including for a letter of numerical
height zero. We fix the initial small coefficients at the vanishing truncation
loss. In the modular counts the remaining coordinates are pairwise coprime by
the anatomic decomposition and primitivity. There is at least one supported
letter since \(3u<1\).
[Fixed-degree cutoff](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#making-the-low-coefficient-short-in-anatomic-boxes) · [Bounded-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

If all three are supported, lock a positive-weight heavy root on either tall
term as a full modulus and release the other two supported terms by (FC2). The
modulus root must cost at least \(d\), and releasing the two chosen tall roots
by (30) gives \(2d\le1-d\), a contradiction.

If exactly two letters are supported, let \(S\) be the remaining letter, of
mass \(s\le u\). If \(S\) has positive numerical height \(h\), fix its part in
degrees at least \(\lfloor h/(2L)\rfloor\). Its numerical weight is positive
since the complementary weight is at most \(sh/2+o(1)\). It is a full modulus
at cost at most \(s\). Release the two supported terms by (FC2), giving
\(d\le s\le u<d\). If \(S\) has zero height, the two supported terms are the
tall terms. Their signed ratio is counted at cost \(o(L)\) by (SR3), again
impossible.

It remains that exactly one letter is supported. Name the unsupported letters
\(U,V\), of masses \(s_U,s_V\le u\). There is a positive-weight high
coordinate of positive normalized mass on each unsupported letter of positive
limiting height. For example, if \(U\) has height \(h>0\), fix \(V\) in full
and consider the high part of \(U\) in degrees at least
\(\lfloor h/(2L)\rfloor\). Were its largest coordinate mass at most
\((d-s_V)/2+o(1)\), group the high coordinates into blocks of mass at most
\((d-s_V)/2+\xi\) for a sufficiently small fixed \(\xi>0\). Start a new block
when the next coordinate would exceed the cap. Consecutive completed blocks
have combined mass above the cap, so there are boundedly many blocks and at
least one has positive limiting weight. Using that block as modulus and
counting the supported term by (FC2) would cost less than \(d\). Thus a
coordinate of mass bounded below exists on a subsequence, and its degree
threshold ensures positive limiting weight. Denote its mass by \(\beta_U\),
and likewise by \(\beta_V\) when \(V\) has positive height.

If both \(U,V\) have positive height, use the selected coordinate of \(U\) as
modulus and fix the complement of the selected coordinate on \(V\). Release
the selected coordinate on \(V\) along with the supported term by (FC2). The
cost is \(\beta_U+s_V-\beta_V\). Exchanging their roles and adding the bounds
gives
\[
2d\le s_U+s_V\le2u,
\]
a contradiction. If one unsupported letter, say \(V\), has numerical height
zero, then \(U\) and the supported letter are the tall terms. Fix all of \(U\)
apart from its selected positive-weight coordinate. Its fixed part costs
\(s_U-\beta_U\), and (SR3) releases the selected coordinate and the full
supported tall term jointly by their signed ratio. Thus \(d\le s_U-\beta_U\).
On the other hand, use the selected coordinate of \(U\) as modulus and fix
\(V\) in full. Equation (FC2) then gives \(d\le s_V+\beta_U\). Adding again
gives \(2d\le2u\), a contradiction.

This excludes every limiting \(d>1/3\). If the asserted estimate failed for
some fixed \(\eta>0\) on budgets tending to zero, the anatomic truncation and
height-shell decomposition with losses \(o(\lambda)\) would give limiting
boxes with \(L\le\lambda+o(\lambda)\) and \(d\ge1/3+\eta\). The estimates
above exclude those boxes. The endpoint with \(d=1/3\) and all three letter
masses equal to \(1/3\) remains open in this argument.


## Fixed-moment residue branches for the unit budget

Fix the coordinates \(O\) of one letter \(C\) at cost
\(X^{f+o(1)}\). Their product gives a modulus with lower numerical
exponent \(w\). For each of the two released letters \(T\), write
\(L_T\) for its full mass. The checker assigns one of the following
energy exponents \(K_T\).

- A list of distinct early indices \(J_T\), with positive integer
  moment orders \(h_{T,i}\) satisfying
  \(\sum_{i\in J_T}1/h_{T,i}\le1\), gives
  \[
  K_T=2L_T-\sum_{i\in J_T}\alpha_{T,i}
  \qquad\text{if}\qquad
  w\ge h_{T,i}\alpha_{T,i}\quad(i\in J_T).
  \]
  Freeze both copies of the complementary coordinates in the
  congruence energy. Apply (6) with \(k=1\), uniformly in the unit
  twist from those coordinates. Each positive-part term in (6)
  vanishes under the displayed guards. The empty list gives the
  trivial exponent \(2L_T\). An aggregated tail may stay among the
  complementary coordinates at twice its mass.
- A null list gives \(K_T=L_T\) under the guard \(w\ge2L_T\) by
  (7). The argument for (7) uses the full root mass, so this case
  also applies with an aggregated tail.

For each fixed choice of the locked coordinates, the complete
residues of the two released letters agree up to their sign in the
primitive equation. Cauchy–Schwarz on the residue distributions
bounds their pair count by the geometric mean of their congruence
energies. Restoring the lock, (9) gives
\[
D\le f+\frac{K_A+K_B}{2}.
\tag{FM}
\]
The verifier's fixed-moment branch uses (FM) when all guards hold,
with one complementary branch for the weak reverse of each guard.
These branches cover the closed projected box. The integer moment
orders in the verifier are bounded by \(128\).

For two exact unit-budget profiles, lock \(O=\{4\}\) on \(a\), so
\(f=a_4\) and \(w=4a_4\). On \(b\) take the diagonal bound, and on \(c\)
use indices \(1,3\) with both moment orders equal to \(2\). The entries
in each profile column are numerators with the indicated common
denominator.

| Quantity | Denominator \(25986\) | Denominator \(239660\) |
| --- | ---: | ---: |
| Candidate \(D\) | \(13523\) | \(126872\) |
| \(a_4\) | \(4627\) | \(46074\) |
| \(w=4a_4\) | \(18508\) | \(184296\) |
| \(L_b\) | \(5244\) | \(41584\) |
| \(L_c\) | \(11694\) | \(112788\) |
| \(c_1\) | \(4627\) | \(49352\) |
| \(c_3\) | \(6988\) | \(63436\) |
| \(\max(2L_b,2c_1,2c_3)\) | \(13976\) | \(126872\) |

All guards hold strictly. Equation (FM) bounds \(D\) by
\(a_4+(L_b+2L_c-c_1-c_3)/2\), giving \(8757/17324\) and
\(6163/11983\), respectively. Both caps are strictly below \(13/25\).

The universal interpolation also admits a tighter rational envelope.
For fixed \(y=L_T>0\), the convex bound in the proof of (8), divided
by \(y\), is bounded at \(w/y=0,1,3/2,2\) by
\(2,293/250,259/250,1\), respectively. The two interior bounds follow
from (7), \(707^2<2\cdot500^2\), and \(433^2<3\cdot250^2\).
Chords between successive vertices, together with the diagonal
bound for \(w\ge2y\), give
\[
K_T\le\max\left\{
2y-\frac{207}{250}w,\quad
\frac{361}{250}y-\frac{34}{125}w,\quad
\frac{143}{125}y-\frac9{125}w,\quad y
\right\}.
\]
The displayed lines have successively increasing slopes and meet at
\(w/y=1,3/2,2\). The zero-mass case follows from the trivial bound.

## Necessary conditions at the balanced \(1/3\) endpoint

The [outside factor graph](ABC_DIFFUSE_FACTOR_GRAPH.md) treats the
case where the summed nonatomic degree mass is \(1/2\). Fixing
the outside factors on any two letters bounds the number of
completions uniformly, so the resulting graph has nearly quadratic
size on the root-box scale.

Consider a limiting sequence of anatomic boxes with \(L\to0\) and
\[
d=l_A=l_B=l_C=1/3.
\]
Take the exceptional coefficients fixed at the vanishing truncation loss. Write \(\beta_{T,e}=L_{T,e}/L\) for a coordinate's normalized root mass and \(eL\beta_{T,e}\) for its numerical weight. The degree set at each truncation is finite. Coordinates within a letter may be taken pairwise coprime, and the three complete terms are pairwise coprime. We use (FC2), (SR3), the bounded-root release (30), and the one-high scalable congruence.

### Residual masses and their numerical weights

For fixed nonnegative integers \(J,M\), discard on each letter all coordinates of degree at most \(J\), then discard the \(M\) largest remaining coordinates in normalized mass, breaking ties by degree. Write \(r_T(J,M)\) and \(v_T(J,M)\) for the limiting normalized mass and limiting numerical weight of the coordinates left over, along a subsequence on which all the indicated limits exist. These quantities are nonincreasing in \(J,M\). Define
\[
\delta_T=\inf_{J,M}r_T(J,M),\qquad
\omega_T=\inf_{J,M}v_T(J,M).
\]
Thus \(\delta_T\) is the mass that remains after removing the bounded-degree part and arbitrarily many macroscopic coordinates. The corresponding residual numerical weight is \(\omega_T\). These definitions concern the given limiting sequence and its chosen subsequence.

If a full modulus from \(T\) of positive limiting weight costs normalized mass \(\mu\), and its fullness threshold can be chosen of order \(1/L\), then
\[
\frac13\le\mu+\delta_U+\delta_V
\tag{D1}
\]
for the other two letters \(U,V\). Indeed, freeze their residual coordinates at finite \(J,M\). The coordinates of degrees at most \(J\) on each free letter combine into one integer of height \(X^{O_J(L)}\). The \(M\) other retained coordinates on that letter are individually of height \(X^{O(L)}\). Apply (FC2) to this bounded collection, using the signed ratio of the free terms. The coprime ratio determines both complete values and the triple. The normalized cost is \(\mu+r_U(J,M)+r_V(J,M)\). Take \(J,M\) large after taking the small-budget limit. The same proof allows a divisor of a coordinate or a block of coordinates as the modulus, since its prime support stays on \(T\).

If the short term has a logarithmic height gap \(\tau>0\) to the two tall terms with \(\tau/L\to\infty\), apply (SR3) after freezing the same residual coordinates on the two tall letters. The relative discrepancy of the two tall terms is \(O(X^{-\tau})\); the retained tuples have height \(X^{O_J(L)}\), and each exact complete ratio determines both tall values. Therefore
\[
\delta_B+\delta_C\ge1/3.
\tag{D2}
\]
This includes every case where the limiting short height is strictly below \(1\).

If \(\omega_T>0\), the modulus in (D1) can have arbitrarily small normalized cost. Fix \(J,M\) with \(M\) large. Every coordinate left over has normalized mass at most \(1/(M+1)+o(1)\). Their numerical weight is at least \(\omega_T-o(1)\). Keeping only degrees at least \(\lfloor \omega_T/(2L)\rfloor\) leaves positive numerical weight, since lower degrees contribute at most \(\omega_T/2+o(1)\). Group the high remainder into blocks of normalized mass at most \(2/(M+1)+o(1)\). Consecutive completed blocks have combined mass exceeding the cap, so there are boundedly many blocks for fixed \(M\). Some block has positive limiting weight. It is an available full modulus. Since \(M\) is arbitrary, (D1) gives
\[
\omega_T>0\quad\Longrightarrow\quad \delta_U+\delta_V\ge1/3 .
\tag{D3}
\]

### Positive-weight coordinates and consequences

For a positive-height letter \(T\), suppose \(t=1/3-\delta_U-\delta_V>0\). The part in degrees at least \(\lfloor h_T/(2L)\rfloor\), where \(h_T>0\) is the limiting height, has positive numerical weight. If its largest coordinate mass stays below \(t\) by a fixed margin, the same block construction gives a positive-weight full modulus of cost below \(t\), contrary to (D1). Thus there is a high coordinate with limiting mass at least \(t\) and positive limiting numerical weight. More generally, *each* positive-weight coordinate of limiting normalized mass has mass at least \(t\), since it can itself serve as modulus. Its degree is sufficiently large for the stated fullness threshold, including when its normalized mass tends to zero.

In particular, \(\delta_T\le 1/3-t=\delta_U+\delta_V\) in this case, because the macroscopic coordinate can be removed among the \(M\) selected ones. If \(t\le0\) the same triangle inequality follows from \(\delta_T\le1/3\). Thus the three residual masses satisfy triangle inequalities whenever all terms have positive limiting height. When the short height tends to zero, the two positive letters give their triangle inequalities and the inequality for the short letter follows from (D2) and \(\delta_A\le1/3\).

As a consequence, at most one letter has \(\delta_T=0\). Two zeros force three by the triangle inequalities. Three zeros contradict (D2) if there is a short height gap. With all positive heights, three zeros would force a macroscopic high coordinate of normalized mass at least \(1/3\) on each term, contradicting (30).

Write \(\Delta=\delta_A+\delta_B+\delta_C\). If all heights are positive and \(\Delta<1/3\), choose one positive-weight coordinate on each term as above. Their total normalized mass is at least
\[
\sum_T(1/3-\delta_U-\delta_V)=1-2\Delta.
\]
Their degrees tend to infinity, so (30) bounds this sum by \(1-d=2/3\), giving \(\Delta\ge1/6\). This lower bound is automatic if \(\Delta\ge1/3\); if one height is zero, (D2) gives the stronger \(\Delta\ge1/3\). In fact (D2) and (D3) give \(\Delta\ge1/3\) whenever there is a short height gap as specified above or any \(\omega_T>0\).

These statements also give finite positive-weight structure if \(\Delta<1/3\). Every positive-height term then has \(\omega_T=0\) by (D3). Each of its positive-weight coordinates costs at least \(1/3-\Delta>0\), so only boundedly many such atoms can occur. More formally, pass to a further diagonal subsequence for the sets of \(M\) largest coordinates outside degree \(J\), so their pairwise coincidences stabilize and each selected mass and weight has a limit. Any distinct selected coordinates with positive limiting weight each have mass at least \(1/3-\Delta\), hence induce finitely many such coordinate sequences on each letter. The sum of their limiting weights is the whole letter's weight, since \(\omega_T=0\) and fixed degrees contribute zero. In this regime the diffuse mass is carried at scaled degrees \(eL\) tending to zero, in the sense that for every fixed \(c>0\), its mass in degrees \(eL\ge c\) is at most \(\omega_T/c=0\).

### Single-scale refinement at zero short height

Say a positive-height letter concentrates at a single normalized degree if there is \(r>0\) such that for each \(\epsilon>0\), the numerical weight outside \(|eL-r|\le\epsilon\) tends to zero. Such a letter can be written as \(VT^q\), where \(q=\lfloor r/L\rfloor\), the root \(T\) is the product of the coordinates in a slowly shrinking interval about \(r\), and the residual rational \(V\) has height \(X^{o(1)}\). Indeed, the weight outside the interval tends to zero, and the sum of the root masses times the discrepancies \(|e-q|\) inside the interval also tends to zero. The height of \(T\) is \(X^{O(L)}\).

If the short term has limiting numerical height zero, neither tall term can concentrate at a single normalized degree in a balanced limiting sequence. Suppose the first tall term does so. Using a positive-weight full modulus from the second tall term, the scalable one-high count recovers the ratio of the first tall term and the numerical short term at cost \(o(L)\) beyond the modulus. The residual multiplier comes from \(V\) divided by the short term and has height \(X^{o(1)}\). Thus every admissible positive-weight modulus from the second tall term costs at least \(1/3\). Its total root mass is \(1/3\). The high-block argument above forces a single high coordinate of limiting mass \(1/3\). Any positive height on the remaining coordinates would give a cheaper block, since their total mass tends to zero. So this coordinate carries the full tall height and has \(eL\to3\). The second tall term therefore also concentrates at a single normalized degree. Reverse the modular argument to force the same full-mass single-coordinate structure on the first. Both tall residual diffuse masses vanish, contradicting (D2).

Thus a zero-short-height balanced limiting sequence needs at least \(1/3\) diffuse mass on the two tall letters together, and numerical weight on both tall terms spread across more than one limiting normalized degree or escaping to arbitrarily large normalized degrees.

The [four-prime graph argument](ABC_WIDE_COFACTOR_CONGRUENCES.md#strict-residual-mass-at-the-balanced-13-endpoint) gives the stronger necessary condition \(\Delta>1/3\) for every balanced limiting box. The [bounded-variable extension](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13) forces every limiting profile at \(d=1/3\) to be balanced and gives \(\Delta\ge1/2\), or \(\Delta\ge2/3\) if the short height tends to zero. These lower bounds hold already for the nonatomic part of the limiting root mass distributions in the scaled degrees.
