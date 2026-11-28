# Small complementary cuts between a simplex and a score winner

Fix \(0<r<1\) and a radius \(R<\infty\). Let \(A\) be an
ordinary planar partition given almost everywhere by its
adjusted Gaussian winning scores at \(r\), with positive
cell masses. Let \(T\) be a centered regular three-cell
simplex whose relative orientation is chosen so that no
adjusted pair contrast of \(A\) vanishes identically on
any corresponding full supporting line of \(T\).
Only finitely many orientations are excluded.

Write \(d=1_A-1_T\) in indicator-vector notation. There are
an integer \(L\ge1\) and constants \(c,m_0>0\) such that if
\(E\) is a measurable cut preserving the masses of both
inputs, meaning \(\int 1_Ed\,d\gamma_2=0\), and the transfers
on \(E\) occur in \(B_R\), then
\[
 \langle 1_Ed,T_r(1_{E^c}d)\rangle
 \ \ge\ c\,m^{\,2-1/(L+1)},\qquad
       0<m=\gamma_2(E\cap\{A\ne T\})\le m_0 .                \tag{CM1}
\]
The theorem applies to co-maximizers at any positive
correlation. Both input partitions having stationary winner
scores suffices. If every zero of a relevant restricted
contrast is simple, one can take \(L=1\).

## Regret at both inputs

Write \(U^A_i=T_r1_{A_i}\), \(U^T_i=T_r1_{T_i}\) and use
prices \(\beta^A_i\) for \(A\). The simplex selects the
winner of its unadjusted scores, by reflection of its
two cells across every pair line. On the changed set
with \(A=i\) and \(T=j\), put
\[
 M(x)=U_i^A(x)-U_j^A(x)-\beta_i^A+\beta_j^A
                      +U_j^T(x)-U_i^T(x).
\]
Both summands are nonnegative. The balanced cut condition
cancels the prices in its integral. Thus with \(h=1_Ed\)
and \(Q_r(h)=\langle h,T_rh\rangle\),
\[
             \langle h,T_r(d-h)\rangle
                       =\int_{E\cap\{A\ne T\}} M\,d\gamma_2-Q_r(h).
                                                               \tag{CM2}
\]
Twice this quantity is the difference between the sum
of the two input stabilities and that of their two
complementary patches, by direct quadratic expansion.
[Complementary patch identity, G12](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

The simplex score \(U_j^T-U_i^T\) has the sign of the signed
distance to the supporting line for \(i,j\), oriented
toward \(T_j\). Moreover on \(B_R\) its absolute value
bounds the absolute distance below by a positive multiple.
In fact \(T_i\) is the reflection of \(T_j\), each in its
corresponding open half-plane. Pair the two reflected
source points in the Mehler convolution. The difference
has the strict expected sign off the line, and its
normal derivative on the line is strictly positive,
by differentiating the paired kernel at zero normal
coordinate. The quotient by signed distance extends
continuously with this positive value on the line.
Its minimum on the compact disk is positive.

## A sublevel estimate for the combined regret

If \(M(x)\le t\) at a changed point of \(B_R\), the distance
of \(x\) from the relevant reference line is at most
\(C_Rt\). Write \(x_\parallel\) for its orthogonal projection
onto that line. The adjusted \(i,j\) contrast \(F_{ij}\)
of \(A\) is uniformly Lipschitz on \(B_R\), so
\[
                         |F_{ij}(x_\parallel)|\le C'_Rt .
\]
The restriction of \(F_{ij}\) to the line is real analytic
and not identically zero. On the projection of the disk
it has finitely many zeros, each of finite order. If it
has any zero use the maximal order on that segment, and
otherwise use order one. On each small neighborhood of
a zero of order \(l\), factor the restriction as the
\(l\)-th power of distance times a function bounded away
from zero, and use a positive lower bound elsewhere.
The length of the sublevel set on the line is at most
\(C t^{1/L}\), where \(L\) is the largest of the orders
among the three pair lines. Integrating across the
transverse strip and adding over changed types gives
\[
      \gamma_2(\{A\ne T,\ x\in B_R,\ M(x)\le t\})
                           \le C t^{\,1+1/L}                  \tag{CM3}
\]
for all sufficiently small positive \(t\). There are
only finitely many excluded orientations because for
each pair the first nonzero homogeneous Taylor term of
\(F_{ij}(x)-F_{ij}(0)\) is a nonzero polynomial, with
finitely many zeros on the unit circle, and an identically
zero line must have one of these directions. A pair
contrast is nonconstant by Gaussian convolution
injectivity and positivity of its two cell masses.

With \(\alpha=1+1/L\), take \(t=(m/(2C))^{1/\alpha}\) in
(CM3) for small enough \(m\). At least half of the changed
mass then has regret \(M>t\), so its integral in (CM2)
is bounded below by \(c m^{1+1/\alpha}\).
On the other hand \(h\) has norm at most \(\sqrt2\)
pointwise and is supported on a set of Gaussian mass
\(m\) in \(B_R\). The Mehler kernel is bounded on
\(B_R\times B_R\), giving \(Q_r(h)\le C_R m^2\).
Since \(1+1/\alpha=2-1/(L+1)<2\), this term is absorbed
for sufficiently small \(m\). This proves (CM1).

For fixed co-maximizers at a right contact, choose the
generic relative orientation in the right contact comparison.
Then (CM1) holds for every balanced cut whose changed part
is in \(B_R\) and whose complementary changed mass is at
least a fixed positive \(\eta\), after decreasing \(c\).
For small \(m\) use (CM1), and for masses bounded away from
zero use the uniform strict comparison for the fixed
input pair.
[Right contact comparison](THREE_CELL_RIGHT_CONTACT_PATCH_PROPAGATION.md#strict-comparison-at-right-contact)
