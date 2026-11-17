# A deformed first image fails quadratic descent

Keep the original Weyl algebra with \([P,X]=[Z,Q]=h\), and put
\[
 \begin{split}
 r&=X(2-3XQ),\\
 D_0&=\frac{1+3XQ}{2}P-3Q^2Z,\\
 R_\varepsilon&=r+\varepsilon(h)D_0,
             \qquad 0\ne\varepsilon(h)\in hK[[h]],
 \end{split}                                                    \tag{1}
\]
over a characteristic-zero field \(K\). The central scalar series
\(\varepsilon(h)\) may have any positive h-order and arbitrarily many
nonzero h-coefficients. These are actual finite global polynomial operators,
and
\[
                         [D_0,R_\varepsilon]=h.                   \tag{2}
\]
The first image is nonpositional, h-regular, and has Long's classical
first image \(r\). Nevertheless no h-regular finite rational-position
operator commuting with \(R_\varepsilon\) has canonical Long's classical
second image. This excludes every global finite second image with that
reduction, in all differential orders. No homogeneous Euler relation
or scalar fifth-root premise is used.

## A rational quadratic chart and its constants

Write \(L=K((h))\), \(F=L(X,Q)\), and work in the full finite
rational differential algebra \(F[P,Z]\) with the original commuting
momenta. Set
\[
                 w=1-3XQ,\qquad u=\frac{w^2}{Q}.
\]
The rational inverse formulas are
\[
 Q=\frac{w^2}{u},\qquad
 X=\frac{u(1-w)}{3w^2},\qquad
 r=\frac{u(1-w^2)}{3w^2}.                                   \tag{3}
\]
In particular \(F=L(u,w)\), while \(F/L(r,u)\) is quadratic:
\[
                 Q=(3r+u)^{-1},\qquad
                 w^2=\frac{u}{3r+u}.                          \tag{4}
\]
The latter rational function is not a square in \(L(r,u)\), since
its valuation at \(3r+u=0\) is odd. The nontrivial deck map is
\[
 \sigma(w)=-w,\qquad \sigma(Q)=Q,\qquad
                         \sigma(X)=\frac2{3Q}-X.              \tag{5}
\]

Let \(\delta_0\) be the vector field of \(D_0\), so
\([D_0,f]=h\delta_0(f)\). Directly,
\[
 \delta_0(r)=1,\qquad \delta_0(u)=0,\qquad
                       \delta_0(w)=-\frac{3w^3}{2u}.           \tag{6}
\]
Thus, on \(L(u,w)\), it is a nonzero rational multiple of
\(\partial_w\) with \(u\) fixed, and
\[
                              \ker_F\delta_0=L(u).             \tag{7}
\]
This is an exact rational constants calculation, not a choice of one
formal branch of (4).

Extend the two commuting derivations \(\partial_r,\partial_u\) of
\(L(r,u)\) uniquely to its finite separable extension \(F\).
Uniqueness makes their extensions commute. The first is \(\delta_0\).
For the second put
\[
 \Pi_u=h\partial_u
       =-\frac{3X^2Q^2}{2w}P-Q^2Z.                           \tag{8}
\]
It has \([\Pi_u,u]=h\), \([\Pi_u,r]=0\), and
\([D_0,\Pi_u]=0\). Their coefficient matrix is invertible over
\(F\); explicitly
\[
 P=2w(D_0-3\Pi_u),\qquad
 Z=-3X^2D_0-\frac{w(1+3XQ)}{Q^2}\Pi_u.                        \tag{9}
\]
All these coefficients are independent of \(h\). Therefore
\[
              F[P,Z]=F[\Pi_u;h\partial_u][D_0;h\partial_r]
                                                               \tag{10}
\]
is an equality of finite rational differential algebras. This change
does not complete either momentum or change finite differential order.
It also preserves h-regularity: reordering introduces only nonnegative
powers of \(h\).

## The exact finite centralizer

In (10),
\[
 \boxed{\quad
   \operatorname{Cent}(R_\varepsilon)
                    =L(u)[\Pi_u;h\partial_u][R_\varepsilon].
       \quad}                                                \tag{11}
\]
The last variable commutes with the displayed coefficient algebra.
To prove this, write a commuting finite operator as
\(S=\sum_{j=0}^N a_jD_0^j\),
\(a_j\in F[\Pi_u;h\partial_u]\). Since \(r\) commutes with that
coefficient algebra, the coefficient of \(D_0^N\) in
\([R_\varepsilon,S]\) is
\[
                           \varepsilon h\partial_r a_N.
\]
It vanishes. Equations (7)--(8) give
\(a_N\in L(u)[\Pi_u;h\partial_u]\). This coefficient commutes with
\(R_\varepsilon\). Subtract
\(a_N\varepsilon^{-N}R_\varepsilon^N\); its leading \(D_0\)-coefficient
is exactly \(a_N\), so the order drops. Finite repetition proves
one inclusion in (11), including the order-zero case. The other
inclusion follows from the exact commutators above. Highest
\(D_0\)-coefficients also prove uniqueness of the polynomial expansion.

This argument allows arbitrary rational position coefficients and no
bound on the finite momentum order. Division by \(\varepsilon\) occurs
in the generic coefficient field, so regularity must still be audited.

## Parameter poles cannot disguise a different residue

Let \(S\) in (11) be h-regular for the rational Gauss valuation.
Expand it uniquely as a finite sum
\[
                 S=\sum_{i,j}c_{ij}(u,h)\Pi_u^iR_\varepsilon^j,
                             \qquad c_{ij}\in L(u).            \tag{12}
\]
Because \(\varepsilon\in hK[[h]]\), one has
\(\overline{R_\varepsilon}=r\) for every allowed series, regardless of
its positive order. Every coefficient \(c_{ij}\) is h-regular. If their least valuation
were negative, divide by that least power and reduce. The nonzero
residue of the normalized sum would be
\[
                    \sum_{i,j}\overline c_{ij}(u)\pi_u^i r^j.
                                                               \tag{13}
\]
Here \(r\) is transcendental over \(K(u)\), and \(\pi_u\) is an
independent momentum symbol over \(K(X,Q)\). These monomials are
linearly independent over \(K(u)\), so (13) cannot vanish. It would
give a negative valuation for \(S\), a contradiction. The same
argument shows equality of the valuation of (12) with the minimum
coefficient valuation.

It follows in particular that
\[
 S\text{ finite and h-regular},\ [R_\varepsilon,S]=0
       \quad\Longrightarrow\quad
                  \overline S\in K(u)[r,\pi_u].                \tag{14}
\]
Such a residue is invariant under (5). The deck map fixes the
derivations \(\partial_r,\partial_u\), by uniqueness of algebraic
extension, so it fixes their cotangent symbols as well. No assumption
of regular coefficients in the generic centralizer subtraction was
silently made.

## Canonical Long's second image does not descend

The canonical four polynomial images are those of
[the exact Long transcription](../dixmier-audit/verify_long_pbw.py),
`long_images()`, in target order \((R,T,D,S)\). In commuting classical
source variables put
\[
 \begin{split}
 B&=3X^2P+2wZ,\qquad \beta=B-9Q^2,\\
 y&=Q-X\beta/3,\\
 s&=y+3X(1+Xy)^2\beta+3Xy^2(4+3Xy).
 \end{split}                                                   \tag{15}
\]
These compact formulas equal that transcription. The polynomial in
(15) has \(\beta\)-degree two, with
\[
                      [\beta^2]s=\frac{X^3(3QX-2)}3.
\]
Equations (8)--(9) give
\(B=-2w\pi_u/Q^2\) classically. Thus the quadratic momentum
coefficient is
\[
 [\pi_u^2]s
   =\frac{4X^3(3QX-2)(3QX-1)^2}{3Q^4}
   =-\frac{4ruX^2}{3Q^3}.                                   \tag{16}
\]
The nontrivial deck map changes it by
\[
 \sigma([\pi_u^2]s)-[\pi_u^2]s
                          =-\frac{16ruw}{27Q^5}\ne0.           \tag{17}
\]
Therefore \(s\notin K(r,u)[\pi_u]\), and in particular it cannot
satisfy (14). This proves the stated nonexistence of a finite regular
second image. Adding a function of \(r\) to the classical second
image would not change this quadratic coefficient.

## Formal gauge control and the remaining scope

There is nevertheless an explicit formal automorphism
\[
 \alpha=\exp\!\left(\operatorname{ad}((\varepsilon/h)D_0^2/2)\right)
                                                               \tag{18}
\]
of the h-completed polynomial Weyl algebra. Since \(\varepsilon/h\in K[[h]]\), this is a finite polynomial
operator inside the adjoint. Its derivation sends the finite algebra
into \(\varepsilon\) times that algebra: every commutator has a factor
\(h\). It raises valuation by at least \(\operatorname{ord}_h\varepsilon\ge1\),
so the exponential converges h-adically coefficientwise in the polynomial
algebra. Since \([D_0^2,r]=2hD_0\) and the next iterated commutator
vanishes,
\[
                       \alpha(r)=R_\varepsilon,\qquad
                       \alpha(D_0)=D_0.                       \tag{19}
\]
This is a formal derivation exponential. It is not being represented
by conjugation with an h-adic unit, nor asserted to be a finite
polynomial automorphism.

For a concrete commuting control, let \(S_b\) be the globally admitted
exact-\(R\) second image in
[the global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md).
It is finite, commutes with \(r\), and has classical image \(s\).
Then \(\alpha(S_b)\) is an h-regular formal polynomial-coefficient
operator commuting with \(R_\varepsilon\), with that same residue.
The obstruction above proves it has no common finite PBW support.
Thus a finite first image and an actual individual conjugate do not
supply a finite commuting Long second image, even though formal
normalization is explicit.

The nonzero hypothesis is essential: for \(\varepsilon=0\), the
undeformed finite \(S_b\) is a commuting second image with that residue.
The positive-order hypothesis supplies \(\overline{R_\varepsilon}=r\)
and the convergent formal normalization. There is no fixed order bound
on \(\varepsilon\) beyond positivity.

This result concerns precisely the full scalar-series family (1), not arbitrary deformed first
images or all changes of polarization. It excludes no unrelated
matrix carrier or nonconstant matrix correction. It makes no
algebraic coefficient-cover claim, no division-ring mate assertion,
and no DC1/DC2 conclusion.

[Exact rational controls](check_deformed_first_image.py) verify the
position field changes, commuting vector fields, Long transcription,
quadratic coefficient and nonzero deck difference; their
[retained receipt](deformed_first_image_controls.json) pins those
finite identities. The all-order centralizer and pole-noncancellation
arguments are written proofs, not consequences of a finite prefix.
