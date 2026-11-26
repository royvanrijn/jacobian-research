# Mixed-degree saddle windows under top-block dominance

Fix distinct nonnegative integer exponent columns \(b_i\) and write
\(d_i=\sum_\nu b_{\nu i}\). Keep only the columns whose coefficients in
\(F_D=\sum_i c_i X^{D b_i}\) are nonzero, and let \(H\) be their block of
largest total degree. Put \(h=|H|\), let \(s\) bound the total number of
columns, and set
\[
 u_i(y)=e^{-\sum_\nu y_\nu}\prod_\nu y_\nu^{b_{\nu i}},
 \qquad
 q_i(t)=c_i(Dt)^{D d_i}u_i(b_i)^D,
 \qquad Q_H(t)=\max_{i\in H}|q_i(t)|.
\]
Use \(0^0=1\) at coordinate faces.

**Theorem.** There are \(\gamma>0\) and \(D_0\), depending only on the
columns, with the following property. Suppose \(D\ge D_0\), and on
a window of \(h\) positive integer orders the inequality
\[
              \max_{j\notin H}\frac{|q_j(t)|}{Q_H(t)}
                    \le e^{-\gamma D}t^{-(s+3)}                 \tag{1}
\]
holds. The moments of \(F_D\) cannot all vanish on that window. An
empty lower block imposes no condition. Uniform versions hold for
actual integer columns \(a_i=Db_i\) at a real scale \(D\) whenever
their normalized entries belong to \(\{0\}\cup[\epsilon,M]\), their
normalized columns are pairwise separated by a fixed positive
distance, and their total degrees and active coefficients define a
top block \(H\). Here \(\epsilon,M\) and the separation bound are
fixed positive constants. The condition on the coefficients is
evaluated in \(\mathbb C\), or after an embedding of their field of
definition into \(\mathbb C\). The existence of \(D_0\) is asserted.

This proves the all-order window of length \(h\) for every fixed
polynomial under sufficiently large common dilation of its exponent
columns. It also applies uniformly to varying coefficients for which
\(\max_{j\notin H}|c_j|/\max_{i\in H}|c_i|\le e^{C D}\), with any
fixed \(C\), if the gap from the largest normalized total degree
to every lower one is bounded below by a fixed positive number.
For arbitrary coefficients it gives windows of length \(h\) in
the tail after the crossover. This tail conclusion assumes
\(D(d_i-d_j)>s+3\) for the top degree \(d_i\) and each lower degree
\(d_j\).

## Reduction to a homogeneous expansion

Use the exact scaled Gamma integral, the disjoint peak boxes \(V_i\),
and the constant \(\gamma\) in
[the separated-support saddle argument](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#exact-integral-and-separation-of-the-peaks).
In particular the bound \(|P_j/P_i|\le e^{-4\gamma}\) holds on the
local complex polydisc at the \(i\)-th peak, where
\(P_i=u_i/u_i(b_i)\). For real normalized columns \(b_i\) the
modulus bound is interpreted as in the uniform version there.
The quotient power \((P_j/P_i)^D\) is holomorphic on this
polydisc because all differences of actual exponents are integers
and any negative exponent can occur only in a coordinate that is
positive at the \(i\)-th peak. Separation of the columns ensures
uniform peak boxes. A compact positive interval for the nonzero
entries controls the derivatives and the active Hessians as in
[the homogeneous mixed-active theorem](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md).

Divide at order \(t\) by \(Q_H(t)^t(tD)^r\), where \(r\) is the
number of variables. This leaves the integral with coefficient
weights \(q_i(t)/Q_H(t)\). The weights on \(H\) are independent
of \(t\); write them as \(\hat q_i\). Call a top peak eligible
when \(|\hat q_i|\ge e^{-\gamma D}\). The noneligible top boxes,
all lower boxes, and the complement of the union of the boxes
have combined integral \(O(e^{-ctD})\). Indeed at a lower box
its own weight is bounded by (1), and every other contribution
contains the small relative peak power. The same comparison
holds in a noneligible top box. Off all boxes the bound follows
from the uniform loss in each \(P_i\) and its common integrable
envelope. These are the bounds used to obtain the homogeneous
remainder, with \(Q_H(t)\) in place of the full maximum.

On an eligible top box write \(\Phi_i^0(y)\) for the phase formed
using just the top block, and \(\Phi_i(y,t)\) for the full phase.
Use the same local logarithm of \(P_i\) for both and the logarithm
of a perturbation of one for the remaining factor. By the
relative peak gap and (1),
\[
 \sup_{y\text{ in the polydisc}}
   \big|D(\Phi_i(y,t)-\Phi_i^0(y))\big|
                    \le C e^{-4\gamma D}t^{-(s+3)}.             \tag{2}
\]
Indeed the lower coefficients relative to \(\hat q_i\) are at
most \(t^{-(s+3)}\), every lower term is multiplied by a relative
peak power of modulus at most \(e^{-4\gamma D}\), and the
top logarithm's argument has modulus at least \(1/2\) for large
\(D\).

Compare the two local integrals on the **same** contour. Use the
active-coordinate translation to the top phase saddle
\((\theta_i,0)\), together with the side contours that verify
Cauchy's theorem. For the top phase the real-part bounds in the
[local saddle estimate](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#uniform-local-saddle-and-endpoint-estimate)
give an absolute integral bound on these contours
\[
       C (tD)^{-\rho_i}
                    \exp(tD\operatorname{Re}\Phi_i^0(\theta_i,0)),
       \qquad
       \rho_i=\tfrac12\#\{\nu:b_{\nu i}>0\}
                         +\#\{\nu:b_{\nu i}=0\}.
\]
The translated contour uses the quadratic and linear losses in
the active and inactive coordinates respectively. The integrals
on the side contours have a fixed exponential loss relative to
the saddle value. By (2) the multiplier
\(\exp(tD(\Phi_i-\Phi_i^0))-1\) on all these contours has
modulus \(O(e^{-4\gamma D}t^{-(s+2)})\). Thus the difference can
be bounded relative to the exact top saddle exponential without
losing an exponential in \(t\).

The normalized moment on the window therefore has the form
\[
 \begin{split}
   M_t&=\sum_{i\ {\rm eligible}} z_i^t A_i(t)+R_t+
                              \sum_{i\ {\rm eligible}} E_i(t),\\
   |R_t|&\le C e^{-ctD},\\
   |E_i(t)|&\le C e^{-4\gamma D}t^{-(s+2)}
                                |z_i|^t(tD)^{-\rho_i},          \tag{3}
 \end{split}
\]
where \(z_i=\hat q_i e^{D\Phi_i^0(\theta_i,0)}\) and \(A_i\)
are precisely the frequencies and analytic amplitudes of the
homogeneous expansion. Only the terms \(E_i(t)\) change when
the foreign degree blocks are included.

## Stability of the homogeneous window argument

Suppose there are vanishing windows of length \(h\) along a
sequence \(D\to\infty\). Pass to subsequences with fixed top and
eligible labels and convergent pure-phase amplitudes. The norm
bound on the \(z_i\) is uniform. Some raw maximizing top label
has \(|z_i|\ge\exp(-C e^{-3\gamma D})\).

For an unbounded sequence of window starts \(m\), use the
relative padding, the maximizing label \(p\), and the filters
\(P(E)(E-z_p)^k\) in
[the unbounded-order argument](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md#unbounded-starting-order).
The degree of each filter is at most \(h-1\). All its coefficients
are uniformly bounded since the \(z_i\) are bounded. Recall
\[
 V_i=\prod_{j\ne i}(|z_i-z_j|+|z_i|/m),\qquad
 W_p=\max_i |z_i|^m(mD)^{-\rho_i}V_i.
\]
The normalizing denominators in the filters have magnitude
comparable to \(W_p(|z_p|/m)^{k-|I|+1}\) where \(I\) is the
selected cluster and \(0\le k<|I|\). In particular they are
bounded below by a constant multiple of \(W_p\), with constants
allowed to depend on the bounded cluster along the subsequence.

If \(|z_i|\ge 1/4\), then \(V_i\ge (4m)^{-(h-1)}\).
At each of the finitely many shifts in the filters, the \(i\)-th
error in (3) is therefore bounded by
\[
       C e^{-4\gamma D}m^{-(s+2)}m^{h-1} W_p=o(W_p).
\]
The finitely many powers from shifting \(|z_i|^m\) are uniformly
bounded. If \(|z_i|<1/4\), compare to the raw maximizing top label,
which gives
\[
 W_p\ge c_0(mD)^{-\max_j\rho_j}
                   m^{-(h-1)}e^{-C m e^{-3\gamma D}} .
\]
Here the missing power of \(mD\) up to a fixed exponent is
absorbed by \(e^{-4\gamma D}\) and the decay \(4^{-m}\).
Again the filtered \(E_i\) are \(o(W_p)\). The remainder \(R_t\)
is exactly of the previously controlled type. Hence the new
errors vanish in the normalized limits, which give the same
impossible positive-amplitude Vandermonde equations.

For a bounded sequence of starts, take the start fixed.
The [bounded-start argument](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md#bounded-starting-order)
orders the affine magnitude exponents of the top saddle terms
in the variable \(D^{-1}\). A strictly winning group owns at
least as many consecutive orders as its size. Its normalized
saddle terms remain bounded and all other groups tend to zero.
By (3), each added error tends to zero after that same
normalization because it has the extra exponential factor
relative to its own term. The limiting consecutive
Vandermonde contradiction is unchanged. These two cases prove
the claimed existence of \(D_0\).

For a fixed polynomial, the weight ratios in (1) are bounded
by \(C_1^D(Dt)^{-D\Delta}\) for some \(C_1\) and fixed positive
smallest top degree gap \(\Delta\), whenever the lower block is
nonempty. This proves the all-order consequence. It also proves
the consequence for coefficient ratios bounded by \(e^{CD}\)
under the stated uniform shape and gap bounds. For the tail at
fixed large \(D\), each lower ratio is a constant depending on
\(D\) and the coefficients times \(t^{-D(d_i-d_j)}\).
The strict degree-gap inequality in the theorem ensures that
each such ratio eventually satisfies (1). At every real scale
used here the original moments still have integer exponents.
For fixed coefficients over a characteristic-zero field, use an
embedding of the finitely generated coefficient field into
\(\mathbb C\).

The scope of this transfer is a top block with a root-polynomial
loss in competing degree weights. Windows at a crossover with
competing total-degree blocks of comparable weight remain open.


## All coefficient vectors up to an exponential starting order

There is also a constant \(\eta>0\) depending only on the exponent
columns such that, for sufficiently large \(D\), every positive
window beginning at \(m\le\exp(\eta D)\) and of length equal to the
actual support size has a nonzero factorial moment. This result has
no coefficient-size hypothesis. It applies uniformly on the same
compact normalized shape sets if the gap between any two distinct
normalized total degrees is bounded below by \(g>0\). Equal total
degrees may occur with any multiplicity.

Fix a real constant \(L>r+s+2\). For each degree block \(G\), write
\(Q_G(t)=\max_{i\in G}|q_i(t)|\), with \(q_i\) as before. As a
function of \(x=\log t\), the value \(\log Q_G(t)\) is an affine
line with slope \(D d_G\). If there are \(b\) active degree blocks,
the upper envelope of their lines has at most \(b-1\) change
points. Call an integer \(t\) clean if one block leads and the
gap from it to every other block exceeds \(L\log(tD)/t\).
A nonclean integer is within distance
\[
                     2 L\log(tD)/(D g)                         \tag{4}
\]
of a change point expressed in the \(t\)-variable. To see this,
follow the runner-up line toward its intersection with the
leading line. The intersection lies within
\(L\log(tD)/(tDg)\) in \(\log t\), and some envelope change occurs
at or before the intersection in that direction. The bound in
\(t\) follows by exponentiation for \(D\) sufficiently large.
At an exact tie the integer itself is a change point.

Choose \(\eta\) small enough that (4) is less than \(1/3\) whenever
\(1\le t\le\exp(\eta D)+s\). Then a change point can account for at
most one nonclean integer, so there are at most \(b-1\) such
integers in the window. The clean integers won by a given block
form a consecutive interval within the window. Indeed, between
two such integers that block remains strictly above every other
line, by convexity of its strict winning interval in \(\log t\).
An intermediate nonclean integer would require a change point
within \(1/3\) and hence inside that interval. If the full window
has length \(s'\), the active support size, then at least one
block \(G\) wins at least \(|G|\) clean orders. Otherwise the sum
of clean orders would be at most \(\sum_G(|G|-1)=s'-b\).

Restrict to a window of \(|G|\) consecutive clean orders won by
that block. Normalize using \(Q_G(t)\) and use the homogeneous
phases for the peaks in \(G\). At every one of these orders,
\[
       \max_{j\notin G}\frac{|q_j(t)|}{Q_G(t)}
                       \le U(t):=(tD)^{-L/t}.
\]
The proof on the translated contour still gives (2), with right
side \(C e^{-3\gamma D}\), since every foreign coefficient has
relative modulus at most \(e^{\gamma D}\) in an eligible local
phase. By reducing \(\eta\), this gives the same bound for
\(E_i(t)\) in (3). In fact \(t e^{-3\gamma D}\) is bounded by
\(C e^{-\gamma D}t^{-(s+2)}\) when \(t\le\exp(\eta D)+s\) and
\((s+3)\eta<2\gamma\), and \(e^{-\gamma D}\) is sufficient in the
stability proof.

A foreign box now contributes a polynomial error. On its box
the integrand before taking the \(t\)-th power has modulus at
most \(U(t)+C e^{-4\gamma D}\). Since \(U(t)\ge e^{-L}D^{-L}\)
for \(t\ge1\), reducing \(\eta<4\gamma\) gives
\[
        (U(t)+C e^{-4\gamma D})^t
                  \le 2U(t)^t=2(tD)^{-L}
\]
for large \(D\). The noneligible boxes within \(G\) and the
complement still have an exponential loss. Thus use the same
expansion (3), with \(C(tD)^{-L}\) added to \(R_t\).
In the unbounded-start filter argument, the maximizing
denominator is bounded below by
\(c(mD)^{-r}m^{-(s-1)}\), because
\(m e^{-3\gamma D}=o(1)\). Since \(L>r+s+2\), the additional
filtered error divided by that bound tends to zero. In the
bounded-start argument a peak with raw weight of modulus one
has order-in-\(D^{-1}\) at most \(r\). The strictly winning
group has no larger order, so the polynomial remainder again
vanishes in the normalization. The homogeneous contradictions
exclude a vanishing window.

These two results leave open the crossover regime for
unrestricted coefficients and starts growing faster than this
exponential bound.

