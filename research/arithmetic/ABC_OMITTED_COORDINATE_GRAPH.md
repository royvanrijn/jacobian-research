# Omitted-coordinate graphs in the \(abc\) decomposition

The powerful-number graph applies to the dyadic boxes in the \(abc\)
exceptional-set problem. The transfer uses the modulus given by one complete
summand and a common power among selected coordinates of the other two
summands. The box comparison below gives a strict gain over the four methods
used in the first version of Browning, Lichtman and Teräväinen and the box
estimates used in the later joint paper with Bernert.
[Powerful-number graph, Proposition 10](../../unpublished/powerful-number-representations/main.tex) · [Browning–Lichtman–Teräväinen, §3](https://arxiv.org/html/2410.12234v1#S3) · [Bernert–Browning–Lichtman–Teräväinen, §§2–4](https://arxiv.org/html/2410.12234v2#S2)

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
