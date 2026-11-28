# Three mismatch types on each reference mirror

Fix \(0<r<1\), let \(A\) be a balanced planar Gaussian
noise-stability maximizer with equal stationary prices, and
write \(U_i=T_r1_{A_i}\). Let \(w_1,w_2,w_3\) be unit vectors
with pairwise angles \(120^\circ\). Assume
\[
             U_i=U_j\quad\text{on }\mathbb R w_k
                    \qquad (\{i,j,k\}=\{1,2,3\}).           \tag{RM1}
\]
The lines are the pair lines of the reference simplex with sector
centers \(w_i\). The hypothesis holds for any maximizing compact
discrepancy from that simplex.

If \(A\) is nonsimplex, then on each of the three lines, every
unordered pair of distinct labels occurs on reflected point pairs
on a set of positive Gaussian measure. The three occurrence
measures for all three mirrors have a common positive lower bound
over such maximizers on a compact positive-correlation band.
Neither assertion bounds the radius of a compact discrepancy.

## A centered reflection with at most two mismatch types

Let \(R\) be reflection in any line through the origin, and
let \(H\) be either open half-plane. Write \(f(x)\) for the
vector of the three indicators and, on \(H\), put
\(b(x)=f(x)+f(Rx)\), \(d(x)=f(x)-f(Rx)\). Summing the four
Mehler interactions between \(x,Rx\) and \(y,Ry\) gives
\[
 S_r(A)=\frac12\int_{H^2}
  \left[(M_r(x,y)+M_r(x,Ry))\,b(x)\cdot b(y)
       +(M_r(x,y)-M_r(x,Ry))\,d(x)\cdot d(y)\right]
           d\gamma_2(x)\,d\gamma_2(y).                      \tag{RM2}
\]
Swapping the two labels of a reflected pair preserves all masses,
leaves \(b(x)\) fixed and reverses the sign of \(d(x)\).
The odd kernel in (RM2) is strictly positive on \(H^2\), since
reflection reverses one normal coordinate in the positive cross
term of the Mehler exponent.

Suppose at most two unordered mismatch types occur, with a mismatch
meaning \(d(x)\ne0\). If there are two types, they share one
label. Orient the two incidence vectors out of that common
label. Their mutual dot product is one, and their squared norms
are two. In the one-type case choose either orientation.
Swapping the pairs of either sign therefore attains the dot-product
maximum for each pair of points in \(H^2\) simultaneously. Strict
positivity implies that a global maximizer must already use one
common sign relative to the chosen orientations. Otherwise the
cross term between its two sign sets strictly improves.

In particular, if the type \(i,j\) occurs and there are at most
two types, the remaining label \(k\) has
\(1_{A_k}(x)-1_{A_k}(Rx)\) of one sign almost everywhere on \(H\).
If the other type is absent, this difference is zero. With one
other type, its points all put the label \(k\) on the same side
by the common-sign conclusion.

## The line equalities determine the tilted transverse first moments

Put \(\sigma^2=1-r^2\) and \(\mu=\gamma_{\sigma^2 I_2}\).
The score equality in (RM1), in coordinates \(s\) tangent to
\(\mathbb Rw_k\) and \(v\) normal to it, gives
\[
     \int_{\mathbb R}(1_{A_i}-1_{A_j})(s,v)
                   e^{-v^2/(2\sigma^2)}\,dv=0
                        \quad\text{for almost every }s.    \tag{RM3}
\]
Indeed the Gaussian convolution of this bounded transverse
average in the tangential variable is zero and its Gaussian
Fourier multiplier has no zeros. For a compact discrepancy
the shared exterior rays yield (RM1) and equal prices by
analytic continuation and decay of the discrepancy.
[The compact-discrepancy slice argument (CD1)](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

Let \(m_i=\int x1_{A_i}\,d\mu\), so \(\sum_i m_i=0\).
Multiplying (RM3) by the tangential coordinate against its
Gaussian density gives \(w_k\cdot(m_i-m_j)=0\).
Since \(w_i-w_j\) is perpendicular to \(w_k\), each
cyclic difference \(m_i-m_j\) is parallel to \(w_i-w_j\).
The sum of the three cyclic differences is zero, and the
first two corresponding edge directions are independent.
Their proportionality factors must therefore agree, giving
\(m_i=cw_i\) for a common scalar \(c\).

At the origin all three adjusted scores tie. Their gradients
are distinct there by the maximal-tie gradient separation
at the Gaussian origin. The conditional derivative identity
gives \(\nabla U_i(0)=r m_i/\sigma^2\), hence \(c\ne0\).
Each pair contrast vanishes on its reference line, with
nonzero gradient in its normal direction. Factoring the
contrast as the signed normal coordinate times its averaged
normal derivative shows that the germs of the three strict
winners are exactly the reference simplex when \(c>0\)
and its antipode when \(c<0\). In particular the type
\(i,j\) occurs at reflection in \(\mathbb R w_k\) on an open
set near the incident germ ray.
[Central maximal-tie gradient separation](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)
[The pair contrast factorization (CG1)](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius)

Now fix one of the reference mirrors and let \(n\) be a unit
normal into \(H\). Since \(m_k\) is tangent to the mirror,
\[
 0=m_k\cdot n=
   \int_H (x\cdot n)\big[1_{A_k}(x)-1_{A_k}(Rx)\big]d\mu(x).
                                                               \tag{RM4}
\]
If this mirror has at most two mismatch types, the difference
has one sign by (RM2) and the occurrence of type \(i,j\).
The weight is strictly positive on \(H\). Thus the spectator
cell \(A_k\) is symmetric about the mirror.

## The symmetric spectator determines the split

Freeze \(A_k\) and let \(F=1_{A_i}-1_{A_j}\).
The split-dependent part of stability is
\(\langle F,T_rF\rangle/2\). With \(a=r/\sigma^2\), the
joint Gaussian density in tangential and normal coordinates
gives
\[
 \langle F,T_rF\rangle
  =\sigma^2\iint e^{a(ss'+vv')} F(s,v)F(s',v')\,d\mu\,d\mu .
\]
By (RM3) the term with \(e^{avv'}\) replaced by one vanishes.
For the remaining term, absorb \(\operatorname{sgn}(v)\) into
each factor \(F\). The kernel becomes
\[
 \sigma^2e^{ass'}(e^{avv'}-1)\operatorname{sgn}(vv'),
\]
which is strictly positive off the mirror in either variable.
Choose the sign of the normal coordinate so that the central
germ has label \(i\) on the positive side of its pair ray.
Replacing \(F\) by \(\operatorname{sgn}(v)|F|\) is balanced
because the spectator is symmetric. It replaces the signed
factors by their absolute values and adds a nonnegative
term for the transverse constant kernel. The latter follows
from positive semidefiniteness of the longitudinal exponential
kernel. All pairings are finite by the Gaussian bounds.
The correctly oriented set has positive measure by the germ,
so any wrongly oriented set gives a strict improvement.
[The same split comparison at a compact discrepancy (CD5)](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#a-spectator-bisector-forces-the-simplex)

Thus an actual maximizer already has the two reflected labels
on opposite sides of the mirror. Since the spectator is
symmetric, reflection exchanges \(i,j\) globally. The
reflection reduction makes the spectator a graph with a
finite analytic height on the entire open half-line.
The exact simplex or antipodal germ supplies a straight
segment of this graph. Analyticity makes the graph a line
on the whole half-line. Its slope and intercept are those
of the germ, so the partition is a simplex. This proves
the stated mismatch restriction.
[Analytic graph for a single reflection (RF1)](THREE_CELL_REFLECTION_REDUCTION.md)

## Uniform occurrence on compact correlation bands

If occurrence masses approach zero, pass to a subsequence
fixing the mirror and the unordered pair. Rotate the reference
lines to fix the \(w_i\), and use strong compactness of balanced
planar maximizers as the correlations converge in the band.
Gaussian smoothing then gives locally uniform score convergence,
so (RM1) passes to the limit and the unadjusted scores still
give its winners. The uniform local simplex comparison keeps
a sequence of nonsimplex maximizers a fixed positive distance
from the simplex orbit. Thus the limit is nonsimplex.
The occurrence measures pass to the limit by strong convergence,
since changing labels on a set of mass \(\eta\) changes the
event of a given unordered reflected pair by at most \(2\eta\).
The vanished type would contradict the restriction just proved.
[Compactness and uniform local comparison](POSITIVE_CORRELATION_CONTINUATION.md)
