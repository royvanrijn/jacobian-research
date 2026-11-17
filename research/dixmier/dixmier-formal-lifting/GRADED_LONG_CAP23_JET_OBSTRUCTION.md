# The even graded Long family has no degree-23 full quantum jet

Work over a characteristic-zero field with the standard Moyal product
on source coordinates \((X,Q,P,Z)\), with
\(\{P,X\}=\{Z,Q\}=1\). Use the canonical global polynomial Long tuple
\[
                         f=(R,T,D,S).
\]
Its exact compact construction is
\(\texttt{dixmier/dixmier-audit/verify_long_pbw.py::long_images()}\).
These are full global polynomials, not merely the pair-only
residue-transport inputs: their degrees are \(3,15,23,11\).
[The graph-curvature verifier](check_long_graph_curvature.py)
independently checks all six Poisson relations using Fraction arithmetic.

**Theorem.** There is no even equivariant full scalar Weyl jet
\[
 F_i=f_i+h^2g_i+h^4k_i+O(h^6),\qquad
             \deg g_i,\deg k_i\leq23,                     \tag{1}
\]
satisfying all six canonical commutators modulo \(h^7\), for the
grading specified below. All four images, including the first
image \(R_h\), may deform. Consequently an even equivariant full
Long lift with a common finite spatial-degree bound must have
that bound at least 24.

This is a complete finite-family obstruction. It does not exclude
non-equivariant or odd-jet degree-23 lifts, or any degree-24 lift.

## The grading and the complete finite coefficient spaces

Assign
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).                    \tag{2}
\]
The compact Long formulas verify homogeneity directly:
\(\alpha=1-3XQ\) has weight zero,
\(\beta=3X^2P+2\alpha Z-9Q^2\) has weight two,
\(y=Q-X\beta/3\) has weight one, and \(Xy\) has weight zero.
Every canonical contraction loses source weight three, so
multiplication by the corresponding power of \(h\) preserves
the grading in the Moyal product.

For an image of weight \(w_i\), the complete monomial basis for
its coefficient at \(h^{2m}\), under spatial cap \(C\), consists of
\[
 X^aQ^bP^cZ^d,\quad
 a=b+4c+2d+6m-w_i\geq0,\qquad a+b+c+d\leq C.              \tag{3}
\]
No coefficient shape or momentum-order restriction is added.
At \(C=23\), the first and second correction counts, in order
\((R,T,D,S)\), are
\[
                (54,80,102,71),\qquad(20,34,47,29).       \tag{4}
\]

More generally (3) forces divisibility by
\(X^{\max(0,6m-w_i)}\). Thus an equivariant series with common
spatial cap \(C\) has \(2m\leq(C+w_i)/3\): bounded spatial
support forces termination in \(h\) in this restricted family.
At the first correction, \(r_2,s_2,t_2,d_2\) are respectively
divisible by \(X^7,X^5,X^4,X^2\).
No averaging argument is asserted to put every bounded lift in
this grading.

## The first correction exists, and its entire ambiguity is certified

Write
\[
 E_{ij}=\Pi^3(f_i,f_j)/24.
\]
The coefficient at \(h^3\) in the required commutators gives
\[
       \{g_i,f_j\}+\{f_i,g_j\}=-E_{ij}.                  \tag{5}
\]
Enumerating exactly (3) produces a rational linear system with
307 columns and 2005 rows. It has rank exactly 305.

The retained certificate supplies a nonzero \(305\)-square
minor modulo the prime \(1000003\), of determinant \(199704\).
All its denominators are units modulo that prime. Hence its
rational rank is at least 305. Two exact independent rational
kernel tuples are
\[
 \begin{split}
 K_1&=(0,0,-3R^2,0),\\
 K_2&=(0,R^4,-4R^3S,0).
 \end{split}                                             \tag{6}
\]
They are the Hamiltonian variations from \(R^3\) and \(R^4S\).
Their displayed supports fit (3) at cap 23. Their independence
is immediate from the nonzero \(T\)-component of \(K_2\) and
its absence from \(K_1\). Exact Fraction substitution verifies
both kernel identities. Therefore the rational rank is at most
305 and the kernel is exactly their span, over every
characteristic-zero scalar extension.

The certificate also retains a 59-coefficient rational particular
solution and verifies its substitution into every row of (5).
Its first correction is zero, and its second-position correction is
\[
                     r_2=0,\qquad
             s_2=\frac{49}{2}X^5-\frac{165}{4}X^6Q.        \tag{7}
\]
The particular \(t_2,d_2\) have respectively 9 and 48 terms,
degrees 11 and 19, and momentum orders one and three.
Every solution of (5) is this particular solution plus
\(\alpha K_1+\beta K_2\). In particular every first correction
\(r_2\) in this complete graded cap-23 family vanishes.
This is not an assumption fixing the first image.

## Removing the two first-jet gauges preserves the required two-jet cap

In the abstract target Weyl algebra, \(R,S\) commute.
For
\[
                         G=\alpha R^3+\beta R^4S
\]
the exact polynomial target momentum shear is
\[
 R\mapsto R,\quad S\mapsto S,\quad
 D\mapsto D-h^2G_R,\quad T\mapsto T+h^2G_S.                \tag{8}
\]
Its inverse has the opposite parameters. The derivatives are
functions of the commuting positions, so this preserves every
Weyl relation. Applied to a jet, it can remove (6), with a choice
of signs.

It is essential to check the cap beyond its first coefficients.
Since \(r_2=0\), the \(h^4\) coefficients of \(R,S,T\) acquire
no contribution from (8). The only extra \(h^4\) coefficient of
\(D\), up to its constant parameter, is
\[
                     R^3s_2+\frac18\Pi^2(R^3,S).          \tag{9}
\]
The second term is required: the target product \(R^3S\) is
evaluated by the source star product. Its order-\(h\) correction
vanishes since \(\{R,S\}=0\), while its order-\(h^2\) correction
is the term in (9). Higher contractions do not affect this
two-jet calculation. Both terms in (9) have degree at most
16, and the replay checks their exact combined support.

The \(\alpha R^2\) and \(\beta R^4\) position powers have no
ordering correction at this order. Thus (8) preserves the
degree-23 cap through \(h^4\) and the grading. Every hypothetical
jet (1) may therefore be reduced to the retained first particular
solution. No assertion that this gauge preserves a complete
degree-23 series at later orders is used.

## The second correction is impossible

With \(g\) now fixed to that particular solution, the coefficient
at \(h^5\) is
\[
 \begin{split}
 \{k_i,f_j\}+\{f_i,k_j\}=-H_{ij},\qquad
 H_{ij}={}&\{g_i,g_j\}\\
 &+\frac{\Pi^3(f_i,g_j)+\Pi^3(g_i,f_j)}{24}
                      +\frac{\Pi^5(f_i,f_j)}{1920}.
 \end{split}                                             \tag{10}
\]
In particular the uncured \(S,T\) defect of the \(h^2\) truncation is
\[
                          H_{ST}=-\frac{49}{4}X^{10}R^2.
                                                               \tag{11}
\]
Equation (10) allows every \(r_4,t_4,d_4,s_4\) in the complete
second basis in (4), including 20 first-image coefficients.
It has 130 columns, 1207 rows and 96 nonzero right-hand-side
monomials.

A retained \(131\)-square minor of its AUGMENTED matrix has
determinant \(331880\) modulo \(1000003\). Hence the augmented
rational rank is at least 131, whereas the coefficient rank is
at most the number of unknowns, 130. The rational system is
inconsistent, and remains inconsistent after every
characteristic-zero scalar extension. This is a nonzero-minor
certificate with a rational rank upper bound, not an inference
from modular inconsistency alone.

The exact gauge reduction above covers all first-jet parameter
choices. No later coefficient of the formal series can change
the coefficient at \(h^5\). This proves the theorem.

## Retained data, replay and scope

[The witness](graded_long_cap23_witness.json) retains the minor row
and column indices, prime and determinants, together with all
59 coefficients of the particular first correction.
[The verifier](check_graded_long_cap23.py) reconstructs both
complete rational systems from the canonical global Long input.
It verifies the particular and both rational kernels exactly,
and checks the two selected determinants with bare integer
modular elimination. FLINT was used only to discover the
witnesses; the retained replay does not require it.
[The receipt](graded_long_cap23_controls.json) records current
input pins and the precise scope.

The result concerns EVEN Moyal symbols in the grading (2).
An odd first jet contributes its own quadratic term to the
next curvature. Non-equivariant first-jet kernel components can
also enter the nonlinear second equations. Neither is removed
here by a bounded normalization. The result supplies no
ungraded or odd-jet exclusion, no cap-24 exclusion, and no
DC1/DC2 resolution. It does not follow from an exact-\(R\)
assumption or a finite-order scalar mate classification.

Independent whole-proof replay, formal verification, external
review and novelty priority remain unclaimed.
