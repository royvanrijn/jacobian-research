# Kernel-gcd refinement of actual square-quotient cycle weights

Let an actual selected simple bipartite graph have squarefree positive
left kernels b<=K1, right kernels d<=K2, and one fixed positive
representation on every edge:
\[
 \kappa_1 a^2b^3+\kappa_2 c^2d^3=n,\qquad
 \gcd(\kappa_1a^2b^3,\kappa_2c^2d^3)=1,\quad a\le R_1,\ c\le R_2.
                                                               \tag{1}
\]
The positive coefficients are fixed while the graph is counted.
Assume n>=2, all caps between1 and n, and
\[
                         2(R_1R_2)^2<n.                       \tag{2}
\]
Only simple six-cycles and selected actual labels are considered.
Any usual labelled-cycle convention changes an absolute constant.

At a shared-right path b1--d--b2 with b1!=b2 put
\[
 \alpha=\gcd(a_1,a_2),\quad a_i=\alpha r_i,\quad
 g=\gcd(b_1,b_2),\quad b_i=g u_i,\quad
 \ell=(c_1^2-c_2^2)/\alpha^2.                                \tag{3}
\]
Below sf(M) is the positive squarefree kernel of a nonzero integer:
the unique squarefree e with |M|=e v^2 for a positive integer v.
Define
\[
 e'=\operatorname{sf}(\ell/g^3),\qquad
 {\cal W}=\alpha^3\sqrt{g/e'}.                               \tag{4}
\]
Square roots in weights are ordinary positive real numbers.

**Theorem 1.** All the quantities in (3)--(4) are defined on every
selected shared-right path, with ell a nonzero integer and g^3|ell.
Uniformly in each fixed positive coefficient pair, for every epsilon>0,
\[
 \sum_{\text{shared-right paths}}{\cal W}
             \ll_\epsilon n^\epsilon R_1^2R_2,               \tag{5}
\]
\[
 \sum_{C_6}\ \sum_{\text{shared-right paths in }C_6}{\cal W}
             \ll_\epsilon n^\epsilon K_1R_1^2R_2.            \tag{6}
\]
Consequently for every real Y>0,
\[
 \#\{C_6:\text{some shared-right path has }{\cal W}\ge Y\}
             \ll_\epsilon n^\epsilon K_1R_1^2R_2/Y.          \tag{7}
\]
There are side-swapped versions with the indices interchanged.
No full six-cycle estimate or improvement of an all-target
prescribed-sum exponent is asserted.

## 1. Integral lift, exact cube divisibility and primitive variables

Actual within-edge coprimality implies that every root, kernel and
coefficient in an edge is a unit modulo n. It also gives
gcd(alpha,kappa2*d)=1 and gcd(g,kappa2*d)=1. Subtracting (1) at the
two edges proves alpha^2 divides c1^2-c2^2 and
\[
 \kappa_1r_1^2b_1^3-\kappa_1r_2^2b_2^3
                              +\kappa_2\ell d^3=0.           \tag{8}
\]
If c1=c2, equality of the positive square-times-squarefree-cube
expressions a1^2*b1^3=a2^2*b2^3 forces b1=b2 by prime-valuation
parity. Thus ell!=0. Since g^3 divides both first terms of (8) and
gcd(g,kappa2*d)=1, it divides ell. After its cancellation,
\[
 \kappa_1r_1^2u_1^3-\kappa_1r_2^2u_2^3
                         +\kappa_2(\ell/g^3)d^3=0.           \tag{9}
\]
The coefficients in (9) are NONZERO SIGNED integers. All three
kernel variables are positive; u_i are squarefree and
gcd(u1,u2)=1. Also each u_i is coprime to d, by actual cross-edge
coprimality. In particular gcd(u1,u2,d)=1, precisely the primitive
condition for the cubic input. No coprimality between alpha and g
is assumed or needed.

We use Heath-Brown, Counting Square-full Solutions to x+y=z,
arXiv2601.07817v1, Theorem1.3: primitive integer zeros of
A x^3+B y^3+C z^3=0 in maximum height H number O_epsilon(H^epsilon),
uniformly in arbitrary nonzero signed integer coefficients.
Here H<=max(K1,K2)<=n, so this is O_epsilon(n^epsilon).
The coefficient set is counted before this primitive cubic input is applied.
The global triple-count theorem, fixed-target line theorem, seventh-power
theorem and all-q transfer are not inputs.

## 2. A coefficient count before the cubic theorem

For dyadic A,G>=1 consider alpha in[A,2A), g in[G,2G), and impose e'<=E
for real E>=1. Roots r_i are at most R1/A. For a fixed integer g,
\[
       0<|\ell/g^3|<L_g=R_2^2/(A^2g^3).
\]
The number of DISTINCT possible nonzero signed integer lifts ell/g^3
with squarefree kernel at most E is bounded by
\[
 2\sum_{\substack{1\le t\le E\\t\ {\rm squarefree}}}
          \left\lfloor\sqrt{L_g/t}\right\rfloor
 \le 4\sqrt{L_gE}.                                         \tag{10}
\]
If L_g<1 there are no lifts. Terms with t>L_g vanish, so (10)
remains valid when E is larger than the actual range.
Squarefree uniqueness is available but not essential for this
upper bound; dropping squarefreeness only enlarges the sum.

Summing over integers G<=g<2G, and allowing also non-squarefree g
only for the upper bound, gives
\[
 \sum_{G\le g<2G}\#\{\ell/g^3\}
                    \ll R_2\sqrt E/(A\sqrt G).
\]
Therefore the total possible ordered coefficient data
(g,r1,r2,ell/g^3) number
\[
                \ll R_1^2R_2\sqrt E/(A^3\sqrt G).            \tag{11}
\]
When R1<A the set is empty; all endpoint floors can only reduce
this count. No extra independent alpha or square-root factor is
summed after the coefficient tuple.

For each coefficient tuple (9), the primitive cubic theorem
gives n^epsilon kernel triples(u1,u2,d). Their original kernels
are b_i=g*u_i. These determine at most ONE selected ordered
actual path, because each of its two edges has a fixed
representation. Thus no value of alpha, the square factor in
ell/g^3, or an individual representation is an additional
multiplicity. The bound (11) times n^epsilon is an ACTUAL path
frequency, rather than merely a bound on completions of one slice.

## 3. Weighted sums and cycle recovery

Partition also e' into dyadic ranges [E0,2E0).
There are O(log(n)^3) triples of dyadic ranges: alpha<=R1<=n,
g<=K1<=n, and 1<=e'<=|ell/g^3|<R2^2<=n^2.
The maximum of (4) on such a range is
\[
                 \ll A^3\sqrt{G/E0}.
\]
Apply (11) with E=2E0. The product is
O_epsilon(n^epsilon R1^2R2), independent of A,G,E0.
Absorbing the logarithms by starting with a smaller epsilon
proves (5).

Under (2), either same-side pair of distinct kernels has only
n^epsilon actual common neighbors. This is the exact general
codegree recovery of the coupled-path theorem, Section2:
the reduced unit ratio(a1*c2)/(a2*c1) has height<=R1R2, and
two different ratios in the same square-root class would have
a nonzero cross-difference divisible by n but smaller than n.
After divisor-many factor allocations the exact original equation
recovers the remaining root square and squarefree kernel.
The hypotheses of that theorem, including actual cross-coprimality
and selected labels, are exactly (1)--(2); it is not an extension
of a balanced-only codegree statement.

A specified shared-right path extends to at most K1*n^epsilon
simple six-cycles: choose the third left kernel and its two common
right neighbors using the two codegrees. This proves (6).
Every cycle counted by (7) contributes at least Y to the positive
sum (6), so Markov gives (7) for EVERY Y>0.

## 4. Relation to the current square-quotient weight

Write gamma=gcd(c1,c2) and e0=|ell|/gamma^2, the actual integer
square quotient in the existing coupled-path theorem. Its weight
is W0=alpha^3/sqrt(e0). This e0 is NOT assumed squarefree.

For comparison let e=sf(ell), |ell|=e*s^2 and delta=gcd(g,e).
Since g is squarefree, toggling the prime-exponent parities gives
\[
              e'=ge/\delta^2,\qquad
              {\cal W}=\alpha^3\delta/\sqrt e.               \tag{12}
\]
Actual edge coprimality gives gcd(alpha,gamma)=1. Since gamma^2
divides c1^2-c2^2=alpha^2*ell, it follows that gamma^2 divides ell
and hence gamma divides s. Consequently
\[
             {\cal W}/W_0=\delta s/\gamma
                     \quad\hbox{is a positive integer}.     \tag{13}
\]
In fact actual cross-edge coprimality also gives gcd(g,gamma)=1.
Write |ell/g^3|=e'*v^2. Then gamma^2 divides ell/g^3; since e'
is squarefree, gamma divides v. Consequently
W*/W0=g^2*(v/gamma) is a positive integer multiple of g^2.
The refined weight is therefore at least g^2 times the old one.
In particular merely
replacing the old actual square quotient by sf(ell) already gives
a valid strengthening; (12) gives the further kernel-gcd gain.

Since W0<=W*, equations(5)--(6) bound the old W0 path and cycle
sums with the SAME prefactors. For every Y>0 a cycle with some
W0>=Y also has W*>=Y, so (7) gives the old tail. This recovers
the [square-quotient cycle theorem](ERDOS_943_COUPLED_PATH_ROOT_GCD.md),
Theorem4, by a separate route from its existing proof.

The original positive roots also imply
\[
 e'\le|\ell/g^3|<R_2^2/(\alpha^2g^3),
 \qquad \alpha^4g^2<R_2{\cal W}.                             \tag{14}
\]
Thus (6)--(7) imply an alpha^4*g^2 weighted bound and tail with
prefactor K1*R1^2*R2^2. They also retain the old alpha^4*gamma
corollary through (13) and the old pointwise inequality.
In particular W*>1/R2 on every actual path. As with the earlier
weight, a threshold merely below1 does NOT control all cycles;
a threshold at most1/R2 does.

## 5. Critical and balanced scopes

In the coefficient-shaped critical rectangle of the current
criterion, the general codegree condition is valid for all sufficiently
large n, uniformly in the coefficients: the shaped caps satisfy
2(R1*R2)^2<=32*n^(44/51)<n. The finitely many smaller targets are
absorbed uniformly in the constant, since a nonempty graph has
kappa_i<=n and all vertices and labels are bounded. The exact resource is
\[
 K_1R_1^2R_2=
 16n^{131/153}/(\kappa_1^{19/9}\kappa_2^{14/9}).
\]
Hence the budget n^(26/51+epsilon) is proved on the sector
\[
 {\cal W}\ge
 n^{53/153}/(\kappa_1^{19/9}\kappa_2^{14/9}),                  \tag{15}
\]
with swapped coefficient exponents for the other orientation.
This follows from an unconditional weighted count, not from the
unproved full critical six-cycle hypothesis.

Equations(12)--(15) are the quantitative comparisons. No actual cycles with strict gain are constructed here.

In the balanced box R1=R2=K1=T, the tail is
O_epsilon(n^epsilon*T^4/Y), meeting T^2*n^epsilon on W*>=T^2.
For a necessary size condition, e'>=1 and (14) give
W*<alpha^(8/3)*R2^(1/3). Thus W*>=Y requires
alpha>Y^(3/8)*R2^(-1/8). At coefficient one in (15), this is
alpha>2^(-1/8)*n^(7/68); in the balanced box at Y=T^2 it is
alpha>T^(5/8). These are necessary conditions only, not constructions.

The unselected residual sector remains open. Constants are uniform
for each fixed positive coefficient pair and are not permission
to sum an arbitrary family of coefficient pairs.
