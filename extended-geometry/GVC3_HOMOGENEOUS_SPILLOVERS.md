# Spillovers from the homogeneous three-variable GVC counterexample

## 1. Setup and the complete positive-phase ladder

Work over a characteristic-zero field. Put

\[
\rho=t^2+xy,\qquad A=\rho+x^2,
\qquad C=y\rho^2-2xt^2\rho-x^3t^2,
\]

\[
P_6=AC^2,\qquad
\Delta=4\partial_x\partial_y+\partial_t^2.
\tag{1.1}
\]

The companion note proves

\[
\Delta^{6m}(P_6^m)=0,
\qquad
\Delta^{6m}(x^2P_6^m)\ne0
\quad(m\ge1).
\tag{1.2}
\]

The same endpoint computation gives a stronger multiplier statement. On
\(\rho=1\), let \(u=x^2\), \(B=1+u\), and

\[
K_m(u)=B^m\int_0^1(1-v^2B^2)^{2m}\,dv.
\]

Writing
\[
J_m(B)=\int_0^B(1-w^2)^{2m}\,dw
\]
gives \(K_m(u)=B^{m-1}J_m(B)\), while
\(J_m(1+u)=J_m(1)+O(u^{2m+1})\). Therefore, for
\(1\le\ell\le m\),

\[
\boxed{
\int_{S^2}x^{2\ell}P_6^m\,d\sigma
=\binom{m-1}{\ell-1}C_m,
\qquad
C_m=\frac{2^{2m}(2m)!}{(4m+1)!!}.
}
\tag{1.3}
\]

The pure coefficient, corresponding to \(\ell=0\), is zero. Thus the
counterexample has an entire positive-phase multiplier ladder, not only the
single detector \(x^2\).

More generally, if
\[
R(u)=\sum_{\ell=1}^d r_\ell u^\ell\ne0,
\]
then
\[
\int_{S^2}R(x^2)P_6^m\,d\sigma
=C_m\sum_{\ell=1}^d r_\ell\binom{m-1}{\ell-1}.
\tag{1.4}
\]
The binomial polynomials in \(m\) are linearly independent. Hence every
nonzero \(R\in u k[u]\) detects the witness for all but finitely many
powers.

## 2. Every power \(\Delta^k\) with \(k\ge6\)

For every integer \(k\ge6\), set

\[
\boxed{P_k=\rho^{k-6}P_6,\qquad \Lambda_k=\Delta^k.}
\tag{2.1}
\]

Then \(P_k\) is homogeneous of degree \(2k\). Since \(\rho=1\) on the
sphere, all angular identities in Section 1 are unchanged. Gaussian radial
separation gives, for every \(m\ge1\),

\[
\boxed{\Lambda_k^m(P_k^m)=\Delta^{km}(P_k^m)=0.}
\tag{2.2}
\]

For \(1\le\ell\le m\), the complete scalar detector is

\[
\boxed{
\begin{aligned}
\Delta^{km+\ell}(x^{2\ell}P_k^m)
={}&2^{(k+2)m+\ell}(km+\ell)!(2m)!\\
&\times\frac{(2km+2\ell+1)!!}{(4m+1)!!}
\binom{m-1}{\ell-1}\ne0.
\end{aligned}}
\tag{2.3}
\]

Indeed, the three-dimensional Gaussian moment is

\[
\mathbb E(x^{2\ell}P_k^m)
=(2km+2\ell+1)!!
\binom{m-1}{\ell-1}C_m,
\tag{2.4}
\]

and Wick contraction multiplies it by
\(2^{km+\ell}(km+\ell)!\). Since applying \(\Delta^\ell\) to
\(\Lambda_k^m(x^{2\ell}P_k^m)\) gives the nonzero scalar (2.3), the mixed
output itself is nonzero.

Consequently the generalized vanishing conjecture fails for every operator
\(\Delta^k\), \(k\ge6\), already in three variables.

### 2.1 Scoped minimum for the internal endpoint lift

The phase variable used above is \(u=x^2\), and the Laurent winding is
\(u^{-1}=x^{-2}\). Each endpoint factor
\[
1-t^2(1+x^2)^2
\]
has exactly one factor of \(x\) in the affine quadric ring. Polynomiality
inside the same three homogeneous variables therefore requires endpoint
contact at least two. Contact \(s\) produces degree \(6s\), hence the first
power in this internal one-profile architecture is \(\Delta^{3s}\) with
\(s=2\), namely \(\Delta^6\). This is not a global minimum-order theorem;
it is an exact minimum for this endpoint-contact homogenization mechanism.

## 3. Exact homogeneous dimension threshold

Call \(\operatorname{HGVC}(n)\) the GVC restricted to homogeneous
constant-coefficient operators in \(n\) variables.

> **Theorem 3.1 — homogeneous GVC dimension classification.**
> \[
> \boxed{
> \operatorname{HGVC}(n)\text{ holds if and only if }n\le2.
> }
> \]
> Equivalently, its first failing dimension is exactly three.

For one or two variables, every homogeneous symbol splits into linear forms
after scalar extension, so the translated split-symbol theorem proves GVC
for arbitrary \(P\). Equations (1.1)--(1.2) disprove it in dimension three,
and identity padding handles every larger dimension.

For unrestricted GVC the current dimension ledger becomes

\[
2\le n_{\mathrm{GVC}}\le3,
\tag{3.1}
\]

because GVC(1) is known and unrestricted GVC(2) remains open.

## 4. Quadratic-rank dichotomy for high powers

Let \(q\) be a quadratic form over an algebraically closed
characteristic-zero field, and consider the homogeneous operator
\(q(\partial)^k\).

> **Theorem 4.1 — exact rank threshold for \(k\ge6\).** For every fixed
> integer \(k\ge6\):
> \[
> \boxed{
> q(\partial)^k\text{ satisfies GVC for every polynomial }P
> \iff \operatorname{rank}(q)\le2.
> }
> \tag{4.1}
> \]

If \(\operatorname{rank}(q)\le2\), then \(q\) factors into at most two
linear forms over the algebraic closure. Hence \(q^k\) is a split symbol,
and the split-symbol theorem applies.

If \(\operatorname{rank}(q)\ge3\), choose a nondegenerate
three-dimensional summand and linear coordinates on it in which
\[
q=4XY+T^2+q_{\mathrm{rest}}.
\]
Use \(P_k\) from Section 2 on those three coordinates and make it independent
of all remaining variables. Every term involving
\(q_{\mathrm{rest}}(\partial)\) kills it, so (2.2)--(2.3) give a
counterexample to \(q(\partial)^k\).

Thus rank three, not merely nonsplitting in some large support, is the exact
quadratic obstruction for every power at least six. The cases
\(q(\partial)^k\) with \(1\le k\le5\) are not settled by this theorem.
In particular, the ordinary second-order Laplacian vanishing conjecture is
not touched.

## 5. Rank-one SIC fails first in three pairs

Introduce dual variables \((\xi_x,\xi_y,\xi_t)\) and put

\[
\boxed{
f_k=(4\xi_x\xi_y+\xi_t^2)^kP_k(x,y,t),
\qquad g_\ell=x^{2\ell}.
}
\tag{5.1}
\]

The coefficient tensor of \(f_k\) has separated rank one and bidegree
\((2k,2k)\). Under the contraction map \(\mathcal E_3\),

\[
\mathcal E_3(f_k^m)=\Delta^{km}(P_k^m)=0,
\tag{5.2}
\]

while, for every fixed \(\ell\ge1\) and every \(m\ge\ell\),

\[
\mathcal E_3(g_\ell f_k^m)\ne0
\tag{5.3}
\]

by (2.3). Hence SIC fails on the rank-one Segre cone in three contraction
pairs, in every even balanced bidegree at least twelve.

In at most two pairs, every homogeneous dual symbol splits after scalar
extension, and the split-symbol theorem proves eventual mixed vanishing for
every rank-one form. Therefore:

> **Corollary 5.1 — separated SIC threshold.** The minimum pair dimension
> in which balanced rank-one SIC can fail is exactly three.

This is separate from unrestricted SIC, whose minimum failing pair dimension
is already two and whose known minimal witness is nonseparable.

## 6. Homogeneous GMC and orthogonal nullcones

In real coordinates
\[
x=X+iY,\qquad y=X-iY,\qquad t=T,
\]
let \(G=(X,Y,T)\) be standard Gaussian. Equations (2.2)--(2.4) give

\[
\boxed{
\mathbb E(P_k(G)^m)=0,
\qquad
\mathbb E(x^{2\ell}P_k(G)^m)\ne0
\quad(m\ge\ell).
}
\tag{6.1}
\]

Thus GMC fails already on even homogeneous polynomials in three real
Gaussian variables. Since GMC holds in dimensions one and two, the first
failing dimension of the homogeneous, and even-homogeneous, Gaussian
moments conjecture is exactly three.

There is also an invariant-theoretic consequence. Regard
\(P_k\in\operatorname{Sym}^{2k}(\mathbb C^N)\) by padding unused variables.
If \(P_k\) lay in the complex orthogonal nullcone, Hilbert--Mumford would
supply a one-parameter subgroup for which every weight of \(P_k\) is
positive. The weights of a fixed multiplier are bounded below, so every
orthogonally invariant functional on \(Q P_k^m\) would vanish for large
\(m\). Equation (6.1) contradicts this for \(Q=x^2\). Hence \(P_k\) is
orthogonally semistable, although every scalar power moment vanishes.

Therefore the Gaussian/spherical power-moment invariants do not cut out the
orthogonal nullcone, even on homogeneous even forms. In two dimensions the
one-sided-support theorem gives the opposite conclusion, so the dimension
threshold for this homogeneous moment--nullcone equality is again three.

## 7. Homogeneous Gaussian-to-sphere transfer

Let \(F,H\) be homogeneous polynomials of degrees \(d,e\), respectively,
and let \(U_N\) be uniform on \(S^{N-1}\). For a standard Gaussian
\(G_N=R_NU_N\), radial independence gives

\[
\mathbb E(H(G_N)F(G_N)^m)
=\mathbb E(R_N^{dm+e})
 \int_{S^{N-1}}HF^m\,d\sigma_N.
\tag{7.1}
\]

For even exponent \(2r\),

\[
\mathbb E(R_N^{2r})
=N(N+2)\cdots(N+2r-2).
\tag{7.2}
\]

If \(F,H\) depend only on the first three coordinates, their Gaussian
moment is independent of the ambient dimension. Applying this to
\(F=P_k\), \(H=x^{2\ell}\), gives, for every \(N\ge3\),

\[
\boxed{
\begin{aligned}
\int_{S^{N-1}}x^{2\ell}P_k^m\,d\sigma_N
={}&\frac{(2km+2\ell+1)!!}
 {N(N+2)\cdots(N+2km+2\ell-2)}\\
&\times\binom{m-1}{\ell-1}C_m\ne0,
\end{aligned}}
\tag{7.3}
\]

whereas
\[
\int_{S^{N-1}}P_k^m\,d\sigma_N=0.
\tag{7.4}
\]

Consequently the kernel of spherical mean is not Mathieu on polynomial
functions of \(S^{N-1}\) for every \(N\ge3\), already inside the even
homogeneous subalgebra. Since \(P_k\) and \(x^{2\ell}\) are even, they
descend through the antipodal quotient and give the same conclusion on
every real projective space \(\mathbb{RP}^{N-1}\), \(N\ge3\).

## 8. Compact groups and homogeneous spaces

The previous transfer immediately produces group-theoretic counterexamples.
Let a compact connected group \(G\) act orthogonally and transitively on a
sphere \(S^{N-1}\), with \(N\ge3\). For a base point \(v\), the orbit map
\[
G\longrightarrow S^{N-1},\qquad g\longmapsto gv
\]
pushes normalized Haar measure to normalized spherical measure. Pulling
back \(P_k\) and \(x^2\) therefore gives finite-type functions satisfying
all pure Haar moments zero and every positive mixed Haar moment nonzero.
Thus the Mathieu conjecture fails for every such group and for the sphere
homogeneous space itself.

In particular this gives direct counterexamples for

\[
\boxed{
SO(N)\ (N\ge3),\qquad
SU(N),U(N)\ (N\ge2),\qquad
Sp(N)\ (N\ge1),
}
\tag{8.1}
\]

using their standard real, complex, or quaternionic first-column spheres.
It also transfers to connected covers such as \(\operatorname{Spin}(N)\),
and to any compact homogeneous space carrying a measure-preserving quotient
onto one of these spheres. Real Stiefel manifolds are immediate examples.
The even witness additionally gives all corresponding unoriented/projective
quotients.

This strengthens the earlier isolated \(SU(2)\) and \(SO(3)\) failures: the
homogeneity is exactly what permits dimension padding through Gaussian
radial separation. Any sufficient Abelian conjectural gate used to imply
the displayed group cases must consequently fail in every displayed rank.

## 9. What does not follow

These spillovers must not be confused with the second-order vanishing
conjecture used in the Hessian-nilpotent/Jacobian programme. Here

\[
\Lambda_k^m=\Delta^{km},
\]

whereas the ordinary Laplacian conjecture assumes
\(\Delta^m(P^m)=0\). The displayed \(P_k\) is not asserted to be Hessian
nilpotent. Therefore this construction does not by itself improve the
ordinary-Laplacian dimension bound, prove or disprove a Hessian-nilpotent
statement, or settle GVC(2).

It does, however, show that any route from arbitrary homogeneous
constant-coefficient operators to the Jacobian/Hessian setting must retain
the second-order operator structure. Homogeneity alone is insufficient:
GVC already fails for powers of one nondegenerate quadratic symbol.

## 10. Reproduction

Run

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
python3 scripts/verify_gvc3_homogeneous_spillovers.py
```

The second checker verifies radial padding for \(k=6,7,8,9\), the first
three multiplier levels, rank-one SIC contractions, Gaussian identities,
and sphere-transfer values in ambient dimensions three through eight. The
all-order results are the coefficient, radial, split-symbol, and
Hilbert--Mumford arguments above, not the bounded loops.

## 11. Literature boundary

The endpoint identity originates in Christopher D. Long's
[*Small Counterexamples to the Gaussian Moments Conjecture*](https://arxiv.org/abs/2607.18186).
The binary positive result used in Sections 3--6 is the repository's
split-symbol theorem and its two-dimensional Gaussian lower-face theorem.
Zhao's classical second-order vanishing conjecture is treated in
[*A Vanishing Conjecture on Differential Operators with Constant
Coefficients*](https://arxiv.org/abs/0704.1691). The group-specific Mathieu
programmes for \(SU(N)\), \(SO(N)\), and \(Sp(N)\) include Kevin Zwart's
[*On the Mathieu Conjecture for SU(N) and SO(N)*](https://arxiv.org/abs/2304.02648)
and
[*An addendum on the Mathieu Conjecture for SU(N), Sp(N) and G2*](https://arxiv.org/abs/2504.01516).
No checked source located before this note states the homogeneous padding
consequences of Sections 7--8.
