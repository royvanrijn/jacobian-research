# Homogeneous saddle windows with unequal active-coordinate counts

Let \(b_1,\ldots,b_s\in\mathbb N^r\) be distinct exponent columns
of one total degree. For all sufficiently large integers \(D\),
\[
                         F_D=\sum_i c_i X^{D b_i}
\]
has a nonzero ordinary factorial moment in every positive
window of length equal to its actual support size, for all
coefficient vectors over any characteristic-zero field. The
positive-coordinate counts of the columns can differ.
There is also a uniform threshold for homogeneous integer
matrices \(a_i\) at a real scale \(D\), when \(b_i=a_i/D\), the
positive entries of the \(b_i\) lie in a fixed \([\epsilon,M]\)
with \(\epsilon>0\), and the \(b_i\) are pairwise separated by a
fixed positive distance. This proof gives existence of the
threshold.

The local analysis is that in
[the homogeneous saddle proof](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#a-uniform-analytic-saddle-amplitude).
Here the local powers of the order can differ. The proof uses
relative frequency distances to handle a selected frequency
tending to zero.

## Local terms and their differences

Use the Gamma integral and monomial peaks of the
[separated-support proof](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md).
Discard inactive coefficient columns and fix their number \(s'\).
If there is one active column the result is immediate.
All total degrees agree so one can remove a common scalar at
each order \(m\). Normalize the constants \(q_i\) of the
homogeneous saddle proof by \(\max_i|q_i|=1\). Eligibility
\(|q_i|\ge e^{-\gamma D}\) is independent of \(m\).
After normalization, the moment at positive order \(t\) is
\[
 M_t=\sum_{i\ {\rm eligible}} z_i^t A_i(t)+R_t,\qquad
 A_i(t)=(tD)^{-\rho_i}B_i(tD),\qquad
 |R_t|\le C e^{-ctD},                                      \tag{1}
\]
where \(z_i=q_i e^{D\Phi_i(\theta_i,0)}\) and
\(\rho_i\) is half the number of positive entries plus the
number of zero entries in \(b_i\). All \(z_i\) are nonzero
and bounded above in modulus, and some eligible \(|z_i|\) is
at least \(\exp(-Ce^{-3\gamma D})\). This follows at a
maximizing \(q_i\) from the small critical value of the
perturbed monomial phase. Neither a common normalization by
\((tD)^\rho\) nor a lower bound on the other frequencies is
used.

Apply the sector estimate of the homogeneous saddle proof
separately at each \(i\), with \(\rho_i\). It gives
\(B_i(\lambda)=K_i+O(|\lambda|^{-1/2})\) on a fixed sector,
where \(K_i\) tends uniformly to its positive pure
amplitude as \(D\) grows. Cauchy estimates give the corresponding
derivative bounds on the error. In particular for each fixed
\(h\ge0\), denoting the forward difference by \(\Delta\),
\[
 \begin{split}
   |\Delta^h A_i(m)|&\le C_h(mD)^{-\rho_i} m^{-h},\\
 \frac{m^h\Delta^h A_i(m)}{(mD)^{-\rho_i}}
        &=(-1)^h(\rho_i)^{\overline h}K_i+o(1)
                                   \quad (m,D\to\infty).
 \end{split}                                               \tag{2}
\]
Indeed derivatives of \((tD)^{-\rho_i}\) give the rising
factorial \((\rho_i)^{\overline h}\); the error and its
derivatives have the additional bound \(O((tD)^{-1/2})\).
An \(h\)-fold difference is the integral of the \(h\)-th
derivative over translates by sums of \(h\) parameters in
\([0,1]\). The sector argument applies also to the uniform
version, as its local holomorphic powers are the actual
integral differences \(a_j-a_i\). Thus all estimates can
be taken on the prescribed compact shape set.

Suppose there is a sequence \(D\to\infty\) of failures.
Pass to subsequences to fix the active and eligible labels,
their zero patterns, and the limits of the bounded positive
pure amplitudes. Use \(m\) for the start of the vanishing
window on this sequence.

## Unbounded starting order

First suppose \(m\to\infty\). For eligible indices put
\[
 d_{ij}=|z_i-z_j|+|z_i|/m,\qquad
 V_i=\prod_{j\ne i}d_{ij},\qquad
 W_i=|z_i|^m(mD)^{-\rho_i}V_i .
\]
Choose \(p\) with \(W_p\) maximal, and keep that label fixed
on a subsequence. The set \(I\) consists of the eligible
indices for which \(m|z_j/z_p-1|\) stays bounded. Pass to a
subsequence so that for all other indices it tends to infinity
and the complex numbers \(w_j=m(z_j/z_p-1)\) converge for
\(j\in I\). Write \(P(z)=\prod_{j\notin I}(z-z_j)\)
over eligible labels, and for \(0\le k<|I|\) test the window
with \(L_k(E)=P(E)(E-z_p)^k\), where \(E\) is forward shift.
Its degree is at most \(s'-1\).

For any polynomial \(L\) the exact Newton difference formula is
\[
 L(E)(z_i^t A_i(t))\big|_{t=m}
       =z_i^m\sum_{h=0}^{\deg L}
            \frac{z_i^h L^{(h)}(z_i)}{h!}\Delta^h A_i(m).
 \tag{3}
\]
Normalize the tested values by the nonzero number
\[
       Z_k=z_p^m(mD)^{-\rho_p}P(z_p)(z_p/m)^k .
\]
Since the distances inside \(I\) are \(O(|z_p|/m)\) and the
outside distances from \(p\) are much larger than \(|z_p|/m\),
\[
                |Z_k|\asymp W_p (|z_p|/m)^{k-|I|+1}.
 \tag{4}
\]
The constants here may depend on the subsequence's bounded
cluster size in relative distance, but stay fixed along it.

For an outside index \(i\), the \(h=0\) term in (3) vanishes.
To bound the other terms use (2) and differentiate the
linear-factor product of \(L_k\). Every nonzero term must
differentiate its factor \(z-z_i\). Charge each removed factor
\(|z_i|/m\) from the difference bound to its padded distance,
including the exact factor \(|z_i|/m\) at the removed self
factor. The bound for the contribution of \(i\) becomes
\[
       C W_i\, (|z_i|/m)\, d_{ip}^{k}
                          \big/\prod_{j\in I}d_{ij}.
\]
For outside \(i\), the factors \(d_{ij}\) for \(j\in I\)
are asymptotic to \(d_{ip}\), because all displacements from
\(p\) within \(I\) are \(O(|z_p|/m)\). Dividing by (4) and
using \(W_i\le W_p\) bounds the relative error by a constant
times
\[
 \frac{|z_i|/m}{d_{ip}}
           \left(\frac{|z_p|/m}{d_{ip}}\right)^{|I|-1-k}.
\]
Both displayed ratios tend to zero. For the first, use
\(|z_i|\le |z_p|+|z_i-z_p|\) and the fact that the outside
distance is much larger than \(|z_p|/m\). This estimate
does not assume that \(|z_p|\) is bounded below.

The remainder (1) remains negligible after applying the
filter and dividing by \(Z_k\). The filter coefficients are
uniformly bounded because all frequencies are bounded. A
raw maximizing \(q_i\) gives
\[
 W_p\ge W_i\ge c_0(mD)^{-\max_j\rho_j}
                    m^{-(s'-1)}e^{-C m e^{-3\gamma D}} .
\]
Here its frequency has modulus at least
\(\exp(-Ce^{-3\gamma D})\), and each padded factor is at
least \(|z_i|/m\). The extra factor in (4) is bounded below
since its exponent is nonpositive and \(|z_p|\) is bounded.
The exponential bound in (1) dominates this polynomial
and small-exponential normalization uniformly.

Inside \(I\) one has \(V_i/V_p\) bounded above and below, since
the frequencies differ by \(O(|z_p|/m)\) and their paddings
are comparable. The power ratios \((z_i/z_p)^m\) converge to
\(e^{w_i}\ne0\). Maximality of \(W_p\) then shows that no
index in \(I\) has \(\rho_i<\rho_p\), while indices with
\(\rho_i>\rho_p\) have vanishing relative amplitude. Write
\(I_0=\{i\in I:\rho_i=\rho_p\}\).
On this cluster \(P(z_i)/P(z_p)\to1\), and for every fixed
\(a>0\),
\[
       (z_p/m)^a P^{(a)}(z_i)/P(z_p)\longrightarrow0 ,
\]
by the factor derivative formula and the outside distances.
Thus only the summands with no derivative landing on \(P\)
survive in the normalized (3). By (2), the vanishing window
imposes the limit equations
\[
   \sum_{i\in I_0}
       \kappa_i e^{w_i}
        \sum_{h=0}^k\binom{k}{h}w_i^{\,k-h}
                    (-1)^h(\rho_p)^{\overline h}=0
                         \qquad(0\le k<|I|),                \tag{5}
\]
where \(\kappa_i>0\) is the pure-amplitude limit. The inner
polynomial is monic of degree \(k\), so these equations
include the ordinary Vandermonde equations on the distinct
\(w_i\). Grouping equal \(w_i\) gives weights
\(e^{w_i}\sum_{j:w_j=w_i}\kappa_j\ne0\), a contradiction.

## Bounded starting order

Otherwise take \(m\) fixed. For eligible \(i\), pass to
limits \(v_i=\lim -\log|z_i|/\log D\) in \([0,\infty]\).
At least one \(v_i=0\) by the raw-maximizing bound. In (1)
the magnitude exponents in the variable \(D^{-1}\) are
the affine functions \(\rho_i+t v_i\), at each positive
order \(t\) of the window. Discard the infinite \(v_i\)
and group identical affine functions. If there are \(g\)
groups, the lower envelope can have ties between distinct
groups at at most \(g-1\) distinct order values. A strict
winning interval for any group is an interval of real
orders, so the integer orders assigned to that group are
consecutive. At least \(s'-(g-1)\) window orders are strictly
assigned. Therefore some group is assigned at least as many
orders as its size, since the sum of all group sizes is at
most \(s'\).

All labels in that group have the same \(\rho_i\). On a
subsequence take one fixed label \(p\) of maximal \(|z_i|\)
there, and pass to limits of \(z_i/z_p\) in the group.
At its strictly assigned orders, divide (1) by
\(z_p^t(tD)^{-\rho_p}\). Other groups vanish because their
magnitude exponents have a positive gap at these orders.
The exponentially small \(R_t\) also vanishes, since \(v_p\)
is finite and the denominator is only polynomially small.
The resulting equations are
\(\sum_i \kappa_i \xi_i^t=0\) for sufficiently many
consecutive orders, on a group with limiting frequencies
\(\xi_i=\lim z_i/z_p\) and weights \(\kappa_i>0\).
At least \(\xi_p=1\); delete the zero frequencies and group
the equal nonzero frequencies. The usual consecutive
Vandermonde system contradicts these positive grouped
weights. This rules out all failure sequences. The
characteristic-zero assertion follows by embedding the
field generated by the coefficients into \(\mathbb C\).

## Density one for arbitrary homogeneous zero patterns

Fix the number of variables, the number of columns, and the
positive-coordinate set of each column. Exclude patterns in
which two columns are forced to be the same pure power of
a single variable. Fill the positive entries of each column
independently with a uniform positive composition of the same
total \(H\). The fraction of arrays satisfying the support-size
window at every order tends to one as \(H\) tends to infinity.

Use the [marginal composition bound](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#density-one-on-homogeneous-supports)
with the number of active coordinates of each individual
column. A column with only one positive entry has that entry
equal to \(H\). For a pair of columns, if either has at least
two positive entries, select a coordinate in that column and
condition on the other column. Its marginal probability
bounds the chance that the sup distance of the two columns
is less than \(gH\). Two columns with one positive entry
each are at distance \(H\) under the pattern hypothesis.
The fraction lost to small positive entries or close pairs is
therefore \(O_{r,s}(\epsilon+g+1/H)\). On the remaining
arrays the uniform theorem applies at \(D=H\), with arbitrary
active-coordinate counts. This density statement allows an
unevaluated threshold with no effectivity assertion here.
