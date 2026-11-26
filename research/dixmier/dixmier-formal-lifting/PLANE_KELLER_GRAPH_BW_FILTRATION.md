# Plane Keller graph modules through two filtration degrees

## Source scope

Kanel-Belov and Fedorov,
[The Generalized DC₁ and JC₂, arXiv:2610.33919v2](https://arxiv.org/html/2610.33919v2),
Theorem 1.4(i), excludes a finite lift of Long's four-variable symplectic
map whose reductions implement the prescribed Frobenius-twisted pullback
on the entire center. It covers every finite matrix size over \(A_2\),
using an integral finitely generated \(\mathbb Z[1/30]\)-model and equality
on a Zariski-dense set of closed points. The corresponding exclusion
was already in [v1](https://arxiv.org/html/2610.33919v1), Theorems 1.4(i)
and 4.2. The formal or bounded-polynomial Long quantization problem still
requires a comparison of its induced arithmetic center map.

The new example in v2 Theorem 1.4(ii) gives a smooth closed exact
Lagrangian affine plane
\[
 (v,z)\longmapsto(v,2z+v^2z^2;\,3v^2z^4,2v^3z^3)
 \quad \text{in }T^*\mathbb A^2.
\]
The theorem now covers dimensions at least two, as compared with at least
three in v1. The displayed cotangent projection is generically of degree
two. The plane Keller graph conjecture is still Conjecture 1.5, with
uniqueness in Proposition 3.11. The local obstruction in v2 Lemma 6.1
uses a residue theorem on a smooth projective curve for scalar orders
two and three over characteristic-zero algebraic function fields.
The v1 Lemma 4.8 treated the rational model with time \(v^2/2\). The
more general order scope used in the repository is recorded in the
[coprime obstruction](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md).

## Actual modules and arithmetic conditions

For a plane Keller map \(f\) over a perfect field \(k\) with
\(p=\operatorname{char}k>\deg f\), v2 Proposition 3.2 gives the graph
splitting module \(M_{f,k}\) as a rank-one projective left generator over
\(A=A_1(k)\). Its right action identifies
\[
 B=B_{f,k}=Z\otimes_{Z',(f^{(1)})^*}A'_1(k)
       \simeq\operatorname{End}_{A}(M_{f,k})^{\mathrm{op}},
 \quad Z=k[X,Y],\quad (X,Y)=(x^p,y^p).                          \tag{1}
\]
The Frobenius twist \(f^{(1)}\) acts on the specified center coordinates
by applying coefficientwise Frobenius to \(f\). The upper degree
condition \(p>d\) for all maps of degree at most \(d\) in Conjecture 1.9
replaces \(p>2d\) in v1.

For a left-ideal realization \(I\) of \(M_{f,k}\), put
\(S=k[\xi,\eta]\), \(g=\gcd(\operatorname{gr}I)\) and
\(J=g^{-1}\operatorname{gr}I\) in the Bernstein filtration. Proposition
5.2 gives \(T=\operatorname{gr}B\subset S\),
\[
 \nu=\dim_k S/J,\qquad \lambda=\dim_k S/T=2\nu,\qquad
 F_t B=\{0\}\cup\{b\ne0:\deg_{X,Y}\operatorname{Nrd}(b)\le t\}.   \tag{2}
\]
The reduced norm uses the central simple algebra over \(\operatorname{Frac}Z\)
and total degree uses \(X,Y\) from (1). The formulas are already in v1
(5)–(8), §7.1, and are retained in v2 (6)–(9), §5.1. Given an exact
filtered left-ideal presentation, one can compute \(J\) by obtaining the
full leading ideal, removing its homogeneous gcd, and counting the
standard monomials. Leading symbols of a generating list may require
additional syzygy reductions to give the full leading ideal.

Theorem 1.8 concerns a map over a finitely generated integral
\(R_0\subset\mathbb C\) and reductions at maps to finite \(k_i\) with
\(p_i\to\infty\), faithful in the sense that \(R_0\hookrightarrow
\prod_{\mathcal U} k_i\) is injective. Bounded \(\nu\) on a set in
\(\mathcal U\) is equivalent to bounded \(\lambda\), to an arithmetic
realization as in Definition 1.7 with the same specialized bimodules,
and to a flat finite type family of dimension-two \(D\)-modules supplying
those specializations. This is the prescribed faithful-family
realization. A lift with constant arithmetic support over a dense open
of its model is the separate formulation in Conjecture 1.5.
Conjecture 1.9 asks for the uniform bound over all finite fields and
maps under the given degree limit. For the implication
\(DC'_1\Rightarrow JC_2\), §5.4 uses a bounded faithful family for each
complex map. The generalized \(DC'_1\) remains an assumption in that
conditional implication.

## Exact two-degree criterion

**Theorem PlaneKellerGraphTwoDegreeLM.** Fix a collection of actual
splitting modules as in (1), with their source characteristic and
degree restrictions. Their BW-invariants are uniformly bounded if and
only if some fixed \(n\ge 2\) satisfies, for every graph in the collection,
\[
 \dim_k F_t B-\dim_k F_{t-1}B=t+1,\qquad t\in\{n,n+1\}.        \tag{3}
\]
The same equivalence holds with the collection restricted to a set of
indices in \(\mathcal U\) in Theorem 1.8.

For any unital graded subalgebra \(T\subset S\) of finite colength
\(\lambda\), if \(T_t\ne S_t\) for \(t>0\), then each unordered pair
of positive integers adding to \(t\) contains a degree \(i\) with
\(T_i\ne S_i\). Indeed products from \(S_i S_{t-i}\) span \(S_t\).
The \(\lfloor t/2\rfloor\) pairs are disjoint as pairs of degrees, with
the middle pair interpreted as a singleton. Including degree \(t\)
itself gives \(\lambda\ge\lfloor t/2\rfloor+1\). Consequently
\(T_t=S_t\) for \(t\ge2\lambda\), with \(T_0=S_0\) by unitality.

Conversely \(T_n=S_n\) and \(T_{n+1}=S_{n+1}\) imply \(T_t=S_t\) for
every \(t\ge n(n-1)\). Indeed write \(t=qn+r\), \(0\le r<n\); then
\(q\ge r\), and \(t=(q-r)n+r(n+1)\) gives the required spanning
products. Thus \(\lambda\le\sum_{i=1}^{n(n-1)-1}(i+1)\). Finally each
increment in (3) has dimension \(\dim_k T_t\), proving the equivalence
using (2).

A proof of the remaining graph bound can therefore give \(t+1\)
independent degree-\(t\) leading symbols of actual right idealizer
elements in two degrees chosen solely from \(d\). For a supplied
exact left ideal \(I\), an element of its quotient division ring
is tested by right multiplication against generators of \(I\).
Equivalently one can work directly with (1) and reduced norms.
The resulting common \(n\) for all graphs, or for one faithful
family for each complex map in the conditional implication, is
still open here.

## Exact control

For \(f=(x,y+g(x))\) with \(g\in k[x]\) and \(p>\deg f\), (3) holds
in every positive degree. The Weyl automorphism
\(x'\mapsto x,\ y'\mapsto y+g(x)+c\) realizes the specified center
map with \(c^p=-g^{(p-1)}\). Such a \(c\) lies in \(k\) because
\(g^{(p-1)}\) is a constant and \(k\) is perfect. Its action follows
from the characteristic-\(p\) Ore identity
\((y+g+c)^p=y^p+(g+c)^p+g^{(p-1)}\) (v2 §7.2, with the
Weyl-algebra form also used in Appendix A). Equation (1) identifies the graph module with the
regular left module carrying this right action, so \(I=A\) gives
\(\nu=\lambda=0\). In the opposite direction v2 Example 5.8 shows
that general rank-one projective modules over even one fixed
finite field have unbounded BW-invariants. The example does not
identify its ideals as splitting modules for Keller graphs.

The [retained v2 PDF](sources/kanel-belov-fedorov-2610.33919v2.pdf)
contains the exact statements and proofs used in (1)–(2) and §5.4.
