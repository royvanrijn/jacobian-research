# Mixed factorial windows with unprotected singleton degree blocks

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer exponent
columns with degrees \(d_i=\sum_\nu b_{\nu i}\) and coordinate
sets \(S_i=\{\nu:b_{\nu i}>0\}\). Suppose that for every \(i\)
whose degree occurs at least twice, and every \(j\) of a
different degree, \(S_j\setminus S_i\ne\varnothing\).

**Theorem.** For sufficiently large integers \(D\), the
ordinary factorial moments of \(\sum_i c_i X^{D b_i}\)
satisfy the actual support-size window at every positive order,
uniformly over the coefficients in characteristic zero. The
same existence assertion for a uniform threshold holds for
integer \(Db_i\) at a real scale \(D\) on compact normalized
shape families with positive lower bounds on nonzero entries,
pairwise center distances, and distinct degree gaps.

A degree block with one column can receive perturbations
from any other block, even when their coordinate sets are
nested. For example the degree-one block \(X,Y\) and the
degree-two singleton \(XY\) satisfy the hypothesis. If a foreign
monomial survives at a singleton saddle, its contribution
is included in the exact saddle value and its first order
derivative in the moment index. The finite-difference filter
needs only one power loss there, since it prescribes just
one zero in that degree group.

The [incomparable-coordinate theorem](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md)
supplies the annular interpolation and the repeated-block saddle
estimates. The new input is a weaker first-difference estimate
for variable singleton amplitudes and its use in the same filter.

## The shifted expansion

Use the peak variables \(q_i(t),P_i,Q,d_*\) and the eligible
labels of the incomparable-coordinate proof, with \(t=m+\ell\)
and \(0\le\ell<s'\), where \(s'\) is the actual number of
nonzero coefficients. A block that becomes a singleton after
discarding zero coefficients can retain the old estimates
or use the estimates below. By the
[exponential-range theorem](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order),
only \(m>e^{\eta D}\) needs consideration. The normalized
integral at offset \(\ell\) divides by
\((tD)^r Q^t e^{Dd_*\ell}\). In this normalization the
exterior and noneligible contribution is \(O(e^{-cmD})\).

For a repeated block every foreign quotient has high-order
vanishing at an inactive coordinate, so the expansions and
errors from the incomparable-coordinate proof apply. Write
its mode as
\[
             u_i^m v_i^\ell \sigma_{G_i}(\ell) A_i(m+\ell),
 \qquad
 \sigma_G(\ell)=
 \exp\!\left(Dd_G[(m+\ell)\log(1+\ell/m)-\ell]\right).
 \tag{1}
\]
Here \(u_i=(q_i(m)/Q)e^{D\varphi_i}\),
\(v_i=u_i e^{D(d_i-d_*)}\), and \(\varphi_i\) is the
order-independent block saddle value. The amplitudes have the
usual bounds
\[
 \Delta^h A_i(m)=O((mD)^{-\rho_i}m^{-h}),\qquad
 m^h\Delta^h A_i(m)/(mD)^{-\rho_i}
             =(-1)^h(\rho_i)^{\overline h}K_i+o(1),
 \tag{2}
\]
with \(K_i\) approaching the positive pure-phase amplitude.
For any prescribed \(L\) the local error is bounded by
\(C_L(mD)^{-L}|u_i|^m |v_i|^\ell(mD)^{-\rho_i}\).

For an eligible singleton use the full phase
\[
 \Phi_i(y,t)=\log P_i(y)+D^{-1}\log\!\left(
   1+\sum_{j\ne i}\frac{q_j(t)}{q_i(t)}(P_j(y)/P_i(y))^D\right).
\]
On the local complex box and on a time disc
\(|t-m|\le a m/D\) with \(a>0\) fixed and small, the logarithmic
perturbation is \(O(e^{-3\gamma D}/D)\), with the same bounds
in local spatial derivatives of each fixed order on a smaller
box. Indeed time varies through factors
\((t/m)^{D(d_j-d_i)}\) bounded on that disc. The eligible
bound at \(m\) and the local peak separation
\(|P_j/P_i|\le e^{-4\gamma}\) then apply. The full saddle
at inactive coordinates zero is holomorphic in \(t\) there,
by the active contraction with a uniform inverse Hessian,
and its exact phase value \(\varphi_i(t)\) satisfies
\(\sup|\varphi_i(t)|=O(e^{-3\gamma D}/D)\).
Let
\[
   \psi_i(t)=tD\varphi_i(t),\qquad
   u_i=(q_i(m)/Q)e^{D\varphi_i(m)},\qquad
   v_i=u_i e^{D(d_i-d_*)+mD\varphi_i'(m)} .
 \tag{3}
\]
The extra frequency factor has logarithm
\(O(D e^{-3\gamma D})=o(1)\). By the relative time disc
and Cauchy's formula,
\[
             \psi_i(m+\ell)-\psi_i(m)-\ell\psi_i'(m)
                             =O(D^2e^{-3\gamma D}/m)
 \tag{4}
\]
for every indicated shift.

The local saddle integral, with its exact phase peeled off,
has analytic amplitude \((tD)^{-\rho_i}B_i(t)\), with
\(B_i(t)=K_i(t)+O((mD)^{-1/2})\) on a smaller relative time
disc. The leading Gaussian-endpoint amplitude \(K_i(t)\)
differs from its positive pure-phase value by
\(O(e^{-3\gamma D}/D)\). These bounds follow by the
[local saddle estimate](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#uniform-local-saddle-and-endpoint-estimate)
and the [complex time sector argument](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#a-uniform-analytic-saddle-amplitude).
In the time disc the saddle displacement is small, the real
quadratic and endpoint losses persist under the small rotation
of \(tD\), and the same translated-contour Taylor estimate
bounds the peeled error. Dependence on \(t\) is holomorphic
since the integral can equally be written on the original box
using the chosen logarithm before deformation.

Include the multiplier which is the exponential of the
difference in (4) in the amplitude. The full singleton mode
is exactly (1), with the new \(u_i,v_i\) and this amplitude.
The derivative of the peeled \(B_i(t)\) near the finitely
many shifts is bounded by
\[
   O\!\left(\frac{D}{m}(mD)^{-1/2}+
               \frac{e^{-3\gamma D}}m\right)=O(1/m)
\]
by Cauchy estimates on the error and on \(K_i\) minus the
fixed pure value. The factor \((tD)^{-\rho_i}\) also varies
by \(O(1/m)\) relative over the shifts. Therefore the only
estimates required for the singleton are
\[
 \begin{split}
 A_i(m)&=(mD)^{-\rho_i}(K_i^*+o(1)),\qquad K_i^*>0,\\
 |\Delta^h A_i(m)|&\le C_h(mD)^{-\rho_i}/m
                                   \quad(1\le h\le s'-1).
 \end{split}                                                     \tag{5}
\]
For the last line, each positive finite difference is a
bounded sum of first differences. Their bound also includes
the multiplier from (4). Only existence of a subsequential
pure-amplitude limit is needed in the first line.

All eligible labels thus satisfy the same annular bounds as
before, with \(R_G=e^{D(d_G-d_*)}\):
\[
          c_0 e^{-\gamma D}R_{G_i}\le |v_i|\le C_0R_{G_i}.
\]
A raw maximizing eligible label has
\(|u_i|\ge e^{-C e^{-3\gamma D}}\), including when it is a
singleton. Hence the normalization and error bounds of the
annular filter are unchanged. The repeated-block errors can
lose any fixed power of \(m\), and a singleton has no
approximation error outside its chosen amplitude.

## Filtering singleton modes

Use the relative padding
\[
 d_{ij}=|v_i-v_j|+|v_i|/m,\quad
 V_i=\prod_{j\ne i}d_{ij},\quad
 W_i=|u_i|^m(mD)^{-\rho_i}V_i,\quad W_p=\max_i W_i.
\]
A limiting relative cluster \(I\) around \(p\), defined by
bounded \(m|v_i/v_p-1|\), lies in one degree group by the
annular separation. Prescribe all outside nodes as roots and
\(k\) extra copies of \(v_p\) for \(0\le k<|I|\), exactly as in
[the annular interpolation](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#interpolation-on-separated-annuli).
For each degree group \(G\) the resulting polynomial \(P_G\)
already contains \(\sigma_G(\ell)\) in its coefficients. Thus
a local mode tests to the Newton difference expression
\[
 u_i^m\sum_{h\ge0}
          \frac{v_i^hP_{G_i}^{(h)}(v_i)}{h!}\Delta^h A_i(m).
 \tag{6}
\]
Normalize by the factor used in the incomparable-coordinate
proof, of modulus comparable to
\(W_p (|v_p|/m)^{k-|I|+1}\).

Suppose an outside label \(i\) is a singleton. Its own group
has just the prescribed simple root \(v_i\). By formula (4)
of the annular interpolation, on the circle of radius
\(|v_i|/4\) about \(v_i\) one has
\[
 |P_{G_i}(z)/(z-v_i)|
       \le C\, d_{ip}^{k}\!\!\prod_{\substack{j\notin I\\j\ne i}}\!d_{ij}.
\]
This uses that all roots in the product lie in other degree
annuli, so the annular monomial and constants have
magnitudes comparable to the cross-group distances.
Consequently for \(h\ge1\),
\[
 |v_i|^h |P_{G_i}^{(h)}(v_i)|
       \le C_h |v_i|\,d_{ip}^{k}
                      \!\!\prod_{\substack{j\notin I\\j\ne i}}\!d_{ij}.
\]
The \(h=0\) summand of (6) vanishes. With (5), the remaining
sum has magnitude at most
\[
 C W_i (|v_i|/m)\,d_{ip}^{k}/\prod_{j\in I}d_{ij}.
\]
This is exactly the outside bound of the old filter. After
normalization it tends to zero by \(W_i\le W_p\) and
\(m|v_i/v_p-1|\to\infty\).

If \(p\) is a singleton, then \(I=\{p\}\) and the filter uses
\(k=0\). There is no same-group root in \(P_{G_p}\).
Its derivatives satisfy
\(|v_p|^h|P_{G_p}^{(h)}(v_p)|\le C_h |P_{G_p}(v_p)|\)
by the same circle bound. Thus (5) makes its normalized
positive-difference terms tend to zero, while the \(h=0\)
term tends to \(K_p^*>0\). All outside terms tend to zero,
including those from repeated blocks by (2) and the old
outside estimate. If \(p\) is in a repeated block, its inside
terms are all covered by (2) and give the original
positive-amplitude Vandermonde contradiction. In both cases
the exterior and repeated-block errors vanish by the
estimates in the original annular proof. Hence the tested
window cannot vanish.

The support-size formulation is that of
[Edo and van den Essen, Conjecture 2.8](https://arxiv.org/abs/1304.3956).
For general characteristic-zero fields one embeds the finitely
generated coefficient field into \(\mathbb C\). The remaining
unrestricted case permits foreign monomials that are active
at saddles in a repeated degree block, where a single
difference loss need not control clustered frequencies.

