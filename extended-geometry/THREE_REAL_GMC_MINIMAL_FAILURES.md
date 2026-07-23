# Minimal Gaussian-moment failures in three real variables

## 1. Scope and the first classified island

Let \(X,Y,T\) be independent standard real Gaussians and put

\[
 Z=\frac{X+iY}{\sqrt2},\qquad W=\frac{X-iY}{\sqrt2}.
\]

Then

\[
 \mathbb E(W^rZ^sT^q)
 =\delta_{r,s}r!\,\mathbb E(T^q).
\]

Long's witness is

\[
 P=W+WZ-T^2-\frac32ZT^2-\frac12Z^2T^2.
\]

The global questions are whether five monomials and degree four are minimal,
which supports allow failure, and whether the minimal locus has finitely many
orbits.  This note proves the first exact local classification and specifies
what an exhaustive global certificate must contain.  It does **not** claim
global five-term or degree-four minimality.

Consider the rank-one quadratic ansatz

\[
 P=W h(Z)+v(Z)T^2,\qquad
 h(Z)=1+aZ,\quad
 v(Z)=b_0+b_1Z+b_2Z^2.                         \tag{1}
\]

The coefficient of \(W\) has been normalized to one.  The chart \(a\ne0\)
is the first chart in which the circular pair and the real Gaussian genuinely
interact.

> **Theorem 1 (rank-one classification).**  On \(a\ne0\), the following are
> equivalent:
>
> 1. \(\mathbb E(P^m)=0\) for every \(m\ge1\);
> 2. \(\mathbb E(P^m)=0\) for \(m=1,2,3\);
> 3.
>    \[
>    b_0=-a,\qquad b_1=-\frac32a^2,\qquad
>    b_2=-\frac12a^3.                              \tag{2}
>    \]
>
> Every such polynomial is a GMC failure, since
> \[
> \mathbb E(ZP^m)=m!a^{m-1}\ne0\qquad(m\ge1).       \tag{3}
> \]
> It has exactly five monomials and total degree four.  Consequently five
> terms and degree four are minimal **inside (1)**.

The exact checker is
[`verify_three_real_gmc_rank_one_classification.py`](../scripts/verify_three_real_gmc_rank_one_classification.py).

## 2. Sparse Pareto minimality

The rank-one result extends to every small support in a Gaussian circular
frame.

> **Theorem 2 (sparse Pareto theorem).**  Let
> \[
> P=\sum_{s\in S}c_sW^{i_s}Z^{j_s}T^{k_s},
> \qquad i_s+j_s+k_s\le4,\qquad |S|\le4.            \tag{4}
> \]
> If \(\mathbb E(P^m)=0\) for every \(m\ge1\), then \(P\) satisfies GMC:
> for every polynomial \(Q\),
> \(\mathbb E(QP^m)=0\) for all sufficiently large \(m\).
> Consequently Long's five-term quartic is minimal for the product partial
> order
> \[
> (\deg P,\#\operatorname{supp}P)\le(4,4)           \tag{5}
> \]
> in circular Gaussian coordinates.

The exhaustive exact checker is
[`verify_three_real_gmc_sparse_pareto.py`](../scripts/verify_three_real_gmc_sparse_pareto.py).
It enumerates all \(59{,}535\) raw supports with at most four terms in
\(\mathcal M_4\), quotients by \(W\leftrightarrow Z\), rejects supports
having a moment with a unique coefficient monomial, and computes the
remaining saturated moment ideals over \(\mathbb Q\).

For four-term supports, \(5{,}718\) active coefficient-torus ideals reach the
exact Gröbner pass.  Moments through order \(20\) leave three supports.
Orders \(21\) through \(24\) remove two finite-cutoff artifacts.  Across all
term counts, the genuine active survivors are the following three families,
plus their \(W/Z\) reflections and coefficient rescalings:

\[
\begin{array}{c|c|c}
\text{support}&\text{coefficient condition}&\text{interpretation}\\ \hline
\{T,Z,W\}&1+2ab=0&\text{isotropic linear form}\\
\{Z,WT,W^3\}&P=Z+aWT-\frac12a^2W^3&
 \text{Gaussian triangular shear of }Z\\
\{ZT,Z^2,WT,W^2\}&
 P=ZT+aZ^2+\frac1{2a^2}WT-\frac1{4a^3}W^2&
 \text{orthogonal product with an isotropic factor}.
\end{array}                                        \tag{6}
\]

The second family is

\[
 P=\exp(-aW D)(Z),\qquad
 D=W\partial_T-T\partial_Z.                         \tag{7}
\]

Here \(D(2WZ+T^2)=0\), \(D(W)=0\), and \(W D\) is locally nilpotent and
Gaussian-skew.  Thus its exponential is a polynomial Gaussian-preserving
automorphism.  It transports the one-sided coordinate \(Z\), so it satisfies
GMC.

The last family factors as

\[
 P=
 \left(Z+\frac1{2a^2}W\right)
 \left(T+aZ-\frac1{2a}W\right).                     \tag{8}
\]

The second factor is isotropic and orthogonal to the first.  In an
orthogonal Witt frame, every power \(P^m\) therefore contains the \(m\)-th
power of one strictly one-sided isotropic coordinate.  A fixed \(Q\) can
balance that weight only for finitely many \(m\).

Supports whose moments vanish termwise through the circuit bound are also
strictly one-sided: their nonzero charges \(i_s-j_s\) all have the same
sign.  They satisfy GMC by the same bounded-weight argument.  These
observations prove the all-order and mixed-moment part of the theorem; the
finite Gröbner cutoff alone would not.

The scope of Theorem 2 is exact but bounded.  It does not exclude a dense
cubic, and it does not exclude a four-term failure of degree greater than
four.  It proves that no polynomial can improve on Long simultaneously in
both displayed complexity coordinates.

## 3. Finite coefficient certificate

Direct Wick contraction gives

\[
\begin{aligned}
\mathbb E(P)&=a+b_0,\\
\mathbb E(P^2)&=
2a^2+2ab_0+3b_0^2+2b_1,\\
\frac13\mathbb E(P^3)&=
2a^3+2a^2b_0+3ab_0^2+4ab_1\\
&\hspace{25mm}+5b_0^3+6b_0b_1+2b_2.
\end{aligned}                                      \tag{9}
\]

Saturating these three equations by \(a\), equivalently adjoining
\(\alpha a-1\), gives the lexicographic Gröbner basis

\[
 a\alpha-1,\qquad
 2b_2+a^3,\qquad
 2b_1+3a^2,\qquad
 b_0+a.                                             \tag{10}
\]

Thus the first three moments already give (2) scheme-theoretically on this
normalized chart, where the coefficient of \(W\) is fixed to one.  Adding
any remaining deletion equation
\(a=0,b_0=0,b_1=0\), or \(b_2=0\) gives the unit ideal after saturation.
The equation \(b_2=0\) is also the degree-at-most-three subchart.  This proves
both minimality assertions within the ansatz without relying on a numerical
search.

## 4. All-order certificate

For arbitrary polynomials \(h,v,A\), expansion according to the number of
\(v(Z)T^2\) factors gives

\[
 \frac1{m!}\mathbb E\!\left(A(Z)P^m\right)
 =
 [z^m]\,
 A(z)h(z)^m
 \left(1-\frac{2zv(z)}{h(z)}\right)^{-1/2}.         \tag{11}
\]

For (2),

\[
 v(z)=-\frac a2(1+az)(2+az)
\]

and hence

\[
 1-\frac{2zv(z)}{h(z)}=(1+az)^2.                   \tag{12}
\]

The formal square root has constant term one, so (11) becomes

\[
 \frac1{m!}\mathbb E\!\left(A(Z)P^m\right)
 =[z^m]A(z)(1+az)^{m-1}.                            \tag{13}
\]

Taking \(A=1\) and \(A=z\) proves the all-order pure and mixed statements.
This is the necessary second half of a symbolic certificate: a finite
moment ideal finds a component, while (12) proves that the component survives
all moments.

## 5. Equivalence and the unique template

The complex orthogonal change

\[
 Z'=aZ,\qquad W'=a^{-1}W
\]

preserves the Gaussian contraction pairing
\(\mathbb E(W^rZ^s)=\delta_{r,s}r!\).  Under this change,

\[
 P_a=a\left(
 W'+W'Z'-T^2-\frac32Z'T^2-\frac12Z'^2T^2
 \right).                                          \tag{14}
\]

Therefore the nonzero family (2) is one orbit under complex orthogonal
changes and nonzero scaling of \(P\).  If equivalence is restricted to the
compact real group \(O(3)\), this normalization is too coarse and the
modulus of \(a\) must be treated separately.  A global classification must
state explicitly which of these two equivalence relations it uses.

## 6. Support circuits

For a circular monomial \(W^iZ^jT^k\), record

\[
 \operatorname{charge}(i,j,k)=i-j,\qquad
 \operatorname{parity}(i,j,k)=k\bmod2.              \tag{15}
\]

A Wick contraction can contribute only when the total charge is zero and
the total \(T\)-parity is even.  Thus every support \(S\) has a relation
semigroup

\[
 R_S=\left\{
 n\in\mathbb N^S:
 \sum_{s\in S}n_s\operatorname{charge}(s)=0,\ 
 \sum_{s\in S}n_s\operatorname{parity}(s)=0\bmod2
 \right\}.                                         \tag{16}
\]

The primitive support circuits are the support-minimal elements of the
Hilbert basis of \(R_S\).  They are a useful first combinatorial filter, but
they do not by themselves classify moment cancellation: factorial Wick
weights and collisions between several relations of the same total length
are essential.  Long's support has charge multiset

\[
 \{1,0,0,-1,-2\};
\]

its two zero-charge atoms cancel at the first moment, while the positive and
negative circuits propagate that cancellation through all orders via (12).

## 7. What an exhaustive minimality certificate must do

Fix a degree \(D\), a term bound \(r\), and the circular-coordinate monomial
set

\[
 \mathcal M_D=\{W^iZ^jT^k:i+j+k\le D\}.             \tag{17}
\]

For each support \(S\subset\mathcal M_D\) with \(|S|\le r\), write
\(P_S=\sum_{s\in S}c_s s\) and form

\[
 I_{S,M}=
 \left\langle
 \mathbb E(P_S),\ldots,\mathbb E(P_S^M)
 \right\rangle:
 \left(\prod_{s\in S}c_s\right)^\infty.             \tag{18}
\]

An exact search has four stages.

1. Quotient supports by \(W\leftrightarrow Z\), \(T\mapsto-T\), and monomial
   permutations induced by the chosen equivalence group.  Use (16) to remove
   supports with no relevant circuit collisions.
2. Compute (18) over several finite fields as a fast rejection filter.
   Recompute every claimed exclusion over \(\mathbb Q\); a modular unit ideal
   is discovery evidence, not the final characteristic-zero certificate.
3. For every surviving component, derive an all-order recurrence, diagonal
   identity, or Gaussian integration-by-parts certificate.  A finite cutoff
   alone never proves that all pure moments vanish.
4. Test mixed moments on each genuine null component.  Components with only
   one-sided torus support satisfy GMC and are not failures; a failure needs
   an explicit \(Q\) and infinitely many certified nonzero mixed moments.

The first global targets are:

\[
\begin{array}{c|c|c}
\text{claim} & D & r\\ \hline
\text{degree minimality} & 3 & \text{all supports}\\
\text{term minimality at Long's degree} & 4 & r\le4\\
\text{minimal-template classification} & 4 & r=5.
\end{array}
\]

There are \(20\) monomials in \(\mathcal M_3\) and \(35\) in
\(\mathcal M_4\), so the raw fixed-frame support counts are finite and small
enough for enumeration.  The hard part is not listing supports; it is
certifying the all-order behavior of the surviving coefficient components
and then passing from fixed circular coordinates to \(O(3,\mathbb C)\)
orbits.

## 8. Current conclusion

The five-term quartic is not an isolated lucky coefficient vector in its
natural rank-one envelope: it is the unique nontrivial template, its first
three moments already cut it out, and its all-order survival is explained by
one forced square identity.  This makes the rank-one chart a model for the
global search:

\[
\boxed{\text{support enumeration}
\ \longrightarrow\ \text{saturated moment ideals}
\ \longrightarrow\ \text{all-order identities}
\ \longrightarrow\ \text{orthogonal orbit classification}.}
\]

The bounded sparse rectangle is now closed: there is no failure of degree at
most four with at most four circular-coordinate monomials.  Global
degree-three exclusion remains open because a cubic may have up to twenty
terms.  Global four-term exclusion without a degree bound also remains open.
