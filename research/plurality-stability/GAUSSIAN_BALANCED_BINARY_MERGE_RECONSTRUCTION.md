# Balanced binary merges of Gaussian partitions

Exact pair intersections reconstruct a common balanced \(k\)-partition from
its binary \(q\)-way merges whenever \(q\ne k/2\). At half size there is an
additional sign at every point. For each positive correlation, this sign
allows a family with the simplex pair intersections and strictly larger mean
binary stability than the simplex merges. One specified triple intersection
removes the additional sign and restores the exact reconstruction.

## 1. Reconstruction from pair intersections

Let \(k\ge 3\), \(1\le q\le k-1\), \(\mathcal S=\binom{[k]}q\), and
\(N=|\mathcal S|\). Let \((E_S)_{S\in\mathcal S}\) be measurable sets on one
probability space \((X,\mu)\). Impose the marginal and pair equations
\[
\mu(E_S)=\frac qk,\qquad
\mu(E_S\cap E_R)=\frac{|S\cap R|}{k}\quad(S,R\in\mathcal S).
\tag{1}
\]
For \(q\ne k/2\), there is a unique balanced measurable partition
\((A_i)_{i=1}^k\), up to null sets, such that \(E_S=\bigcup_{i\in S}A_i\)
up to null sets for every \(S\).

Write \(g=(1_{E_S})_S\), \(v_i=(1_{\{i\in S\}})_S\), and let \(V\) have
columns \(v_i\). Equations (1) say
\[
\mathbb E g=\frac1k\sum_i v_i,\qquad
\mathbb E gg^\top=\frac1k VV^\top.
\]
For every vector in \(\ker V^\top\), its inner product with \(g\) has zero
second moment. Testing a basis shows that \(g\in\operatorname{im}V\) almost
everywhere. The same moment identities show
\(\sum_S g_S=Nq/k\) almost everywhere, since that sum has variance zero.
The matrix \(V\) is injective. If all \(q\)-sums of \(a\in\mathbb R^k\) are
zero, comparing two subsets related by swapping \(i,j\) first gives
\(a_i=a_j\), and then the common value is zero.

Here is the Boolean classification for this incidence space. When \(q=1\),
the fixed sum of the \(g_S\) is one, so \(g=v_i\) for the unique selected
coordinate. When \(q=k-1\) the fixed sum is \(N-1\), so the same conclusion
holds for the unique omitted coordinate \(S=[k]\setminus\{i\}\). For
\(2\le q\le k-2\), put \(g=Va\). Comparing two subsets related by a swap
shows \(a_i-a_j\in\{0,\pm1\}\). If the coordinates of \(a\) are all equal
then \(g\) is constant, which contradicts its fixed sum. Otherwise \(a\)
has two values at distance one. If the higher value occurs \(m\) times,
the difference between the highest and lowest attainable \(q\)-sums is
\(\min(q,k-q,m,k-m)\). This difference must be one. Hence \(m=1\) or
\(m=k-1\), giving \(g=v_i\) or \(g=\mathbf1-v_i\), respectively. The
second pattern has sum \(N(1-q/k)\) and passes the fixed-sum test exactly
when \(q=k/2\). This is the degree-one Boolean classification on a Johnson
slice used here in its elementary incidence form. See Theorem 1.2 and
Section 3 of [Filmus and Ihringer, *Boolean degree 1 functions on some
classical association schemes*](https://arxiv.org/abs/1801.06034).

When the only patterns are \(v_i\), take \(A_i=\{g=v_i\}\). The first
moment of \(g\) and the injectivity of \(V\) give \(\mu(A_i)=1/k\). The
patterns are distinct, so the same construction gives uniqueness.

## 2. The half-size signed family and its repair

Suppose \(k\) is even and \(q=k/2\). Let \(P=I-\mathbf1\mathbf1^\top/k\)
on \(\mathbb R^k\), set \(H=\mathbf1^\perp\), and write \(w_i=Pe_i\).
Represent each pattern by a unique \(F=\epsilon w_i\), with
\(\epsilon\in\{+1,-1\}\), so that
\[
g_S=\frac12+\langle F,1_S\rangle.
\tag{2}
\]
The vectors \(P1_S\) span \(H\) by the swap argument. It follows from
(2) that (1) is equivalent to
\[
\mathbb E F=0,\qquad \mathbb E FF^\top=P/k.
\tag{3}
\]
Indeed these moments hold when \(F=w_i\) with probability \(1/k\)
for every \(i\). Taking the first and pair moments of (2) and using the
spanning property gives the converse.

Let \(p_i^\pm=\mu(F=\pm w_i)\). The first equation in (3) says that
\(p_i^+-p_i^-\) is constant in \(i\). Testing the second equation on
\(e_i-e_j\in H\) gives
\((p_i^++p_i^-)+(p_j^++p_j^-)=2/k\) for every \(i\ne j\).
Consequently each sum is \(1/k\). In terms of the negative-pattern
mass \(t=\sum_i p_i^-\), all feasible pattern laws have precisely
\[
p_i^+=\frac{1-t}{k},\qquad p_i^-=\frac{t}{k}
\quad(0\le t\le1).
\tag{4}
\]
Conversely these laws give (3), and thus (1), for any measurable spatial
placement of the signed patterns.

Fix two distinct labels \(a,b\). Choose \(S_1,S_2\in\mathcal S\) related
by swapping \(a\) and \(b\), and choose \(S_3\in\mathcal S\) containing
both. This is possible because \(k\ge4\) and \(q\ge2\). Set
\(\chi_S(i)=2\,1_{\{i\in S\}}-1\). The product
\(\chi_{S_1}\chi_{S_2}\) is \(+1\) outside \(\{a,b\}\) and \(-1\) on
\(\{a,b\}\). Since \(\sum_i\chi_{S_3}(i)=0\), one has
\[
\sum_i \chi_{S_1}(i)\chi_{S_2}(i)\chi_{S_3}(i)=-4.
\]
Equations (2) and (4) now give
\[
\mathbb E\prod_{j=1}^3(g_{S_j}-1/2)=-\frac{1-2t}{2k}.
\]
With the first two moments already fixed, the single equation
\[
\mu(E_{S_1}\cap E_{S_2}\cap E_{S_3})
=\frac{|S_1\cap S_2\cap S_3|}{k}
\tag{5}
\]
sets this centered moment equal to its value at \(t=0\). Thus (5) forces
\(t=0\), and the proof of reconstruction applies.

## 3. The Gaussian objective

For \(0<\rho<1\), let \(X,Y\) be standard Gaussians in \(\mathbb R^d\)
with cross covariance \(\rho I\). Write \(T_\rho\) for the corresponding
conditional expectation operator. For a partition \(A\), let
\(\mathrm{Stab}_\rho(A)=\Pr(A(X)=A(Y))\), with \(A(x)\) denoting
the cell label. The mean binary objective is
\[
B_\rho(E)=\frac1N\sum_{S\in\mathcal S}
\Pr(X\in E_S,\ Y\in E_S).
\]
For genuine merges this is exactly
\[
B_\rho(E)=\frac{q(q-1)}{k(k-1)}
+\frac{q(k-q)}{k(k-1)}\,\mathrm{Stab}_\rho(A).
\tag{6}
\]
Indeed, for a uniform \(S\), the probability it includes both labels
equals \(q/k\) when the labels agree and \(q(q-1)/(k(k-1))\) when they
differ. Thus maximizing \(B_\rho\) with (1), and with (5) at half size,
is equivalent to the original balanced partition optimization in the
same dimension. Equality cases correspond in both directions.

For a half-size family satisfying only (1), equations (2) and (3) give
\[
B_\rho(E)=\frac14+\frac{k}{4(k-1)}
\langle F,T_\rho F\rangle_{L^2(\gamma_d;H)}.
\tag{7}
\]
The coefficient follows from
\[
\frac1N\sum_{S\in\mathcal S}(P1_S)(P1_S)^\top
=\frac{q(k-q)}{k(k-1)}P.
\]
For example, the left side commutes with all coordinate permutations and
annihilates \(\mathbf1\). Its trace is \(q-q^2/k\), fixing the scalar
multiple of \(P\).

## 4. Strict improvement under pair constraints at half size

Let \(k\ge4\) be even, \(d\ge k-1\), and fix \(0<\rho<1\). Take the
standard simplex partition \(T=(T_i)_i\) on \(H\) with
\(T_i=\{x\in H:x_i=\max_\ell x_\ell\}\) up to null sets. Work first on
\(H\) with its standard Gaussian, and put
\[
F_T(x)=w_{T(x)},\qquad U_i(x)=T_\rho 1_{T_i}(x).
\]
For any distinct \(i,j\), there is an open ball \(D\) in the interior of
\(T_i\) with \(U_i+U_j<2/k\) throughout and with all \(k!\) coordinate
permutation images disjoint. To see this, along \(x=-Rw_j\) the
conditional Gaussian mean \(\rho x\) has its \(j\)-th coordinate below
each other coordinate by \(\rho R\). Therefore \(U_j(x)\to0\) as
\(R\to\infty\). The other \(U_\ell(x)\) are equal by permutation symmetry
and tend to \(1/(k-1)<2/k\). Choose a finite \(R\) at which the strict
sum bound holds. By continuity it persists in a neighborhood, which
contains points in the interior of \(T_i\) with all coordinates distinct.
A sufficiently small ball around one such point has the required
disjoint images.

On each \(\pi D\), with \(\pi\) a coordinate permutation, change
\(F_T=w_{\pi(i)}\) to \(F=-w_{\pi(j)}\). Set \(F=F_T\) elsewhere. With
\(t=k!\gamma_H(D)\), each positive pattern loses mass \(t/k\) and each
negative pattern gains mass \(t/k\). Hence (4) and all of (1) hold
exactly. For \(\delta=F-F_T\), the pointwise first-variation integrand
on \(D\) is
\[
\langle -w_j-w_i,T_\rho F_T(x)\rangle
=\frac2k-U_i(x)-U_j(x)>0.
\]
The same holds on every permuted ball by equivariance. Consequently
\[
\begin{aligned}
\langle F,T_\rho F\rangle-\langle F_T,T_\rho F_T\rangle
&=2\langle \delta,T_\rho F_T\rangle+
  \langle\delta,T_\rho\delta\rangle\\
&>0,
\end{aligned}
\]
since \(\langle\delta,T_\rho\delta\rangle
=\|T_{\sqrt\rho}\delta\|_2^2\). Equation (7) proves strict improvement
in mean binary stability. The balls can have arbitrarily small positive
measure, so this is a local obstruction in measure as well. Taking
cylinders in the extra coordinates gives the construction in every
\(d\ge k-1\).

The overlapping pair equations alone therefore cannot prove the simplex
bound for half-size merges. At other sizes, or under the repaired
half-size constraints, the equivalent problem from (6) is a joint binary
stability inequality.

## 5. A uniform improvement on the independent scalar bound

Let \(a=q/k\) and \(J_\rho(a,a)\) be the stability of a Gaussian
halfspace of mass \(a\). For fixed \(k,q,\rho\) there is a positive
dimension-independent gap below \(J_\rho(a,a)\) for all families satisfying
(1) when \(q\ne k/2\), and for half-size families satisfying (1) and (5).
There is also a gap for half-size families satisfying (1) alone when
\(k\ge6\). At \(k=4,\ q=2\) the pair-only bound is attained in every
dimension at least three.

Write \(\delta_S=J_\rho(a,a)-\Pr(X,Y\in E_S)\). The measurable-set form of
the dimension-free robust Borell inequality gives a halfspace \(H_S\)
of mass \(a\) satisfying
\[
d_S:=\gamma_d(E_S\mathbin{\triangle}H_S)
\le C(\rho,a)\,\delta_S^{\alpha_\rho},\qquad
\alpha_\rho=\frac{(1-\rho)(1-\rho^2)}{4(1+3\rho)}>0.
\tag{8}
\]
This is Theorem 1.4 of [Mossel and Neeman, *Robust optimality of Gaussian
noise stability*](https://arxiv.org/abs/1210.4126), in the set formulation
stated as Theorem 1.3 of [Heilman, *A Variational Proof of Robust Gaussian
Noise Stability*](https://arxiv.org/abs/2108.04950). Each \(\delta_S\)
is nonnegative by scalar Borell. Thus any lower bound
\(\max_S d_S\ge e>0\) uniform over halfspaces of the prescribed mass gives
\[
B_\rho(E)\le J_\rho(a,a)
-\frac1N\left(\frac{e}{C(\rho,a)}\right)^{1/\alpha_\rho}.
\tag{9}
\]
The following estimates use only \(k,q\).

**Nonhalf size.** Write \(H_S=\{\langle u_S,x\rangle\ge \tau\}\)
with unit normals and \(\tau=\Phi^{-1}(1-a)\). The score vector has a
Gaussian law determined by its normal Gram matrix, hence also realized
by unit normals in \(\mathbb R^N\). Let \(W\sim\gamma_N\) in that
realization, put \(\eta=1/(4N)\), and set \(R=\max(1,|\tau|/\eta)\).
The event
\(\{|W|\le2\sqrt N,\ \min_S|\langle u_S,W\rangle|>\eta\}\)
has probability at least \(1/2\) by Markov's inequality for \(|W|^2\)
and the union bound for \(N\) standard normal coordinates in
\([-\eta,\eta]\). Dilating this event by \(R\) gives
\[
\Pr\bigl(\min_S|\langle u_S,W\rangle|>|\tau|\bigr)
\ge \tfrac12 R^N e^{-2N(R^2-1)}.
\]
The factor follows from the Gaussian change of variables on the ball
\(|W|\le2\sqrt N\). On the displayed event, the \(H_S\) membership
patterns at \(W\) and \(-W\) are complementary. As \(q\ne k/2\),
at least one has a coordinate sum different from \(Nq/k\). By Gaussian
symmetry and the fixed sum for the \(E_S\), the event has probability
at most \(2\sum_S d_S\). Thus one can use
\(e=R^N e^{-2N(R^2-1)}/(4N)\) in (9).

**Repaired half size.** Every halfspace of mass \(1/2\) is centered.
The product of the three centered indicators associated to \(H_{S_j}\)
from (5) has mean zero, since it changes sign under \(x\mapsto-x\).
The product for \(E_{S_j}\) has mean \(-1/(2k)\). Replacing one centered
indicator at a time changes the mean by at most \(d_{S_j}/4\), since
the two other factors have absolute value \(1/2\). Thus one can take
\(e=2/(3k)\).

**Pair-only half size with \(k\ge6\).** The normal Gram entries of
centered halfspaces satisfy
\[
\langle u_S,u_R\rangle
=\sin\!\bigl(2\pi(\gamma_d(H_S\cap H_R)-1/4)\bigr).
\]
This follows from the planar angular measure of the intersection, and
also holds for repeated or antipodal normals by continuity. Take four
subsets formed by adjoining one of \(a_1,a_2\) and one of \(b_1,b_2\) to
the same fixed set of size \(q-2\), with all four labels distinct and
outside the fixed set. Order these subsets around the square in which
an edge swaps one label. Put \(c=\cos(\pi/q)\). Under the target pair
moments in (1), the four-by-four proposed Gram matrix has edge entries
\(c\) and opposite entries \(\cos(2\pi/q)\). Its quadratic form on
\(z=(1,-1,1,-1)/2\) equals \(2c(c-1)<0\). The true Gram matrix is
positive semidefinite. Each off-diagonal entry differs from its target
by at most \(2\pi(d_S+d_R)\) by the displayed formula and (1).
Since \(\sum_{i\ne j}|z_i z_j|=3\), it follows that
\(12\pi\max_S d_S\ge 2c(1-c)\). This gives positive
\(e=c(1-c)/(6\pi)\). The four-point Gram obstruction is independent of
the ambient dimension.

**Pair-only half size with \(k=4\).** Assign the three complementary
pairs of two-element subsets to three independent centered coordinate
halfspaces and their complements. Every noncomplementary pair intersects
in measure \(1/4\), and each complementary pair has zero intersection.
Thus (1) holds and every term reaches the scalar bound.

For the reconstructed problem, (6) converts (9) into a strict improvement
on the partition upper bound from independent applications of scalar
Borell. The estimate measures an unavoidable joint deficit and does not
identify the simplex value as the sharp bound.
