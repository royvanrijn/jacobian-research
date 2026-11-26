# Homogeneous factorial windows under square-root separation

Let \(a_1,\ldots,a_s\) be distinct ordinary exponent columns of the
same total degree in \(r\) variables. Set \(u_i=(\sqrt{a_{\nu i}})_\nu\),
\(r_i=\#\{\nu:a_{\nu i}>0\}\), and \(C_i=\prod_{a_{\nu i}>0}2u_{\nu i}\).

**Theorem.** For every \(r,s\) there is \(T\) such that the actual
support-size window holds at every positive order if all positive
entries are at least \(T\), all \(\|u_i-u_j\|_2\ge T\) for \(i\ne j\),
and \(C_i\ge T C_j\) whenever \(r_i>r_j\). The threshold is uniform
in the coefficients, over every characteristic-zero field.
The proof establishes existence of \(T\).

Rows constant across columns can first be removed, since their factors
are common moment scalars. For equal active-coordinate counts the
prefactor condition is automatic, without a restriction on row
scale ratios. For example \((H^4+b_iH^3,H^4-b_iH^3,H)\) for distinct
fixed integers \(b_i\) qualifies as \(H\to\infty\). Division by the
largest exponent makes the columns coalesce, whereas square-root
distances grow proportionally to \(H\). For differing active counts
the prefactor condition holds eventually when all positive entries
are between \(A\) and \(KA\) for fixed \(K\) and \(A\to\infty\).

The [homogeneous saddle theorem](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md)
uses compact normalized exponent shapes. In square-root coordinates
the saddle widths stay bounded even when those shapes coalesce. A
shift polynomial applied to the integrand controls the part outside
the saddle boxes even if a selected frequency tends to zero.

## Normalization and local saddles

Work over \(\mathbb C\), delete zero coefficients, and assume a sequence
of failures with \(T\to\infty\), fixed zero pattern and active support
size \(s'\ge2\). Write \(d=|a_i|\) at the current profile. Set
\[
 p_i(y)=\prod_\nu e^{a_{\nu i}-y_\nu^2}(y_\nu/u_{\nu i})^{2a_{\nu i}},
 \quad q_i=c_i\prod_\nu e^{-a_{\nu i}}a_{\nu i}^{a_{\nu i}},
 \quad w(y)=2^r\prod_\nu y_\nu,\quad f(y)=\sum_i q_i p_i(y).
\]
Power factors with zero exponent are one. Scaling \(X_\nu=t y_\nu^2\)
makes the order-\(t\) moment \(t^{r+td}\) times
\(M_t=\int f(y)^t w(y)\,dy\) on the nonnegative orthant. Put
\(\ell_{ij}=\|u_i-u_j\|_2\), \(L=\min_{i\ne j}\ell_{ij}\) and
\(\rho_i=r_i/2+(r-r_i)\). The pure log phase in an active coordinate
has second derivative \(-2-2a_{\nu i}/y_\nu^2\le-2\), with value
and first derivative zero at its peak. In an inactive coordinate
it is \(-y_\nu^2\). Thus
\[
  p_i(y)\le e^{-\|y-u_i\|^2},\qquad
  \int p_i^{t/2} w\,dy\le K C_i t^{-\rho_i}\quad (t\ge1).    \tag{1}
\]
The second estimate integrates the Gaussian bound. An active Jacobian
relative to \(2u_{\nu i}\) is at most \(1+|y_\nu-u_{\nu i}|\), while
the inactive factor \(2y_\nu\) gives \(t^{-1}\). Also
\(|\log(C_i/C_j)|\le2\sqrt r\,\ell_{ij}\), since the function
\(\log(2u)\) on \([1,\infty)\), extended by zero at zero, is
\(2\)-Lipschitz on that domain.

Fix \(R>1\) with \(s e^{-R^2/2}<1\). Let \(V_i\) have active
intervals \(u_{\nu i}\pm R\) and inactive intervals \([0,R]\).
These boxes are disjoint for large \(L\) and the active centers exceed
\(4R\). Use analytic coordinates \(y_\nu\) in active positions and
\(Y_\nu=y_\nu^2\) at inactive endpoints. On polydiscs with radii \(2R\)
and \(2R^2\) respectively, \(p_j/p_i\) is holomorphic: negative powers
occur only in active positions. Its modulus is its value at the real
coordinate moduli \(|y_\nu|\) and \(\sqrt{|Y_\nu|}\), as the exponential
factors cancel. These moduli stay within a fixed distance of \(u_i\),
where \(p_i\ge e^{-K_R}\) by the bounded local pure Hessian and the
linear inactive phase. Thus \(|p_j/p_i|\le e^{-\ell_{ij}^2/2}\) for
large \(L\), by (1).

Label \(i\) eligible if \(|q_i|\ge |q_j|e^{-\ell_{ij}^2/8}\)
for all \(j\). On its box set
\[
 \Phi_i=\log p_i+\log\!\left(1+\sum_{j\ne i}(q_j/q_i)p_j/p_i\right).
\]
The perturbation and its fixed-order derivatives on smaller polydiscs
are \(O(e^{-cL^2})\). The pure active Hessian at the center is
\(-4\,\mathrm{Id}\), and the pure inactive derivative is \(-1\).
A contraction in the active coordinates at inactive coordinates
zero gives the saddle, displaced by \(O(e^{-cL^2})\).
Let \(\varphi_i\) be its value and \(z_i=q_i e^{\varphi_i}\);
then \(\varphi_i=O(e^{-cL^2})\). For noneligible labels take
the auxiliary node \(z_i=q_i\).

For eligible \(i\) the local integral equals
\(z_i^t C_i t^{-\rho_i}B_i(t)\), with
\[
 B_i(t)=K_i+O(|t|^{-1/2}),\qquad
                   K_i=(\pi/2)^{r_i/2}+O(e^{-cL^2}).        \tag{2}
\]
The peeled amplitude is holomorphic in a fixed sector about the
positive axis and the estimates are uniform there for large \(|t|\).
In fact, the active Jacobian divided by \(C_i\) is
\(\prod_{\nu\ {\rm active}}y_\nu/u_{\nu i}\), and the inactive measure
is \(dY_\nu\). Translate only the active intervals to the saddle.
A side contour has a fixed pure phase loss. On the translated box,
the saddle-subtracted phase, also after rotation by a small
\(\arg t\), has real part at most \(-c_1(\|x\|^2+\sum Y_\nu)\) by
the definite active Hessian and negative inactive slopes. The
Taylor error relative to the quadratic-linear phase is
\(O(\|x\|^3+\|x\|\sum Y_\nu+(\sum Y_\nu)^2)\). Integrating the
difference of exponentials with the real-part bound, then using
Gaussian and linear scaling, gives \(O(|t|^{-\rho_i-1/2})\);
the varying Jacobian obeys the same error. This applies the
[homogeneous saddle calculation](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#a-uniform-analytic-saddle-amplitude)
with uniform active Hessian and normalized Jacobian. Relative-disc
Cauchy bounds and integration of derivatives over unit cubes give,
for fixed \(h\) and forward difference \(\Delta\),
\[
 m^{\rho_i+h}\Delta^h(t^{-\rho_i}B_i(t))\big|_{t=m}
       =K_i(-1)^h(\rho_i)^{\overline h}+o(1)\quad(m\to\infty) . \tag{3}
\]

## Unbounded starting orders

Suppose the starting orders \(m\) of the failure sequence tend
to infinity. Set
\[
 d_{ij}=|z_i-z_j|+|z_i|/m,\qquad V'_i=\prod_{j\ne i}d_{ij},
      \qquad W_i=|z_i|^m C_i m^{-\rho_i}V'_i .
\]
Choose a maximum at \(p\). A noneligible \(i\) has a witness \(j\)
with \(|q_j/q_i|>e^{\ell_{ij}^2/8}\), so \(|z_j|>|z_i|\).
Then \(d_{il}\le(2m+1)d_{jl}\) for \(l\ne i,j\), and
\(d_{ij}\le d_{ji}\). The exponent \(m\ell_{ij}^2/8\) dominates
the loss \(2\sqrt r\,\ell_{ij}\) from the \(C\) ratio,
\(O_{r,s}(\log(2m))\) from polynomial factors, and
\(m O(e^{-cL^2})\) from \(\varphi_j\). Hence
\[
                 W_i\le e^{-c_2 m L^2} W_p
                              \quad (i\text{ noneligible}). \tag{4}
\]
In particular \(p\) is eligible.

Pass to a subsequence fixing \(p\) and a cluster \(I\) such
that \(m(z_i/z_p-1)\to v_i\) is finite for \(i\in I\) and
the modulus tends to infinity outside. Put
\(P(z)=\prod_{j\notin I}(z-z_j)\),
\(F_k(z)=P(z)(z-z_p)^k\), \(0\le k<|I|\).
The filter degree is at most \(s'-1\). Normalize \(F_k(E)M_m\)
by \(Z_k=z_p^m C_p m^{-\rho_p}P(z_p)(z_p/m)^k\), where \(E\)
is the forward shift. Its modulus is comparable along the
subsequence to \(W_p(|z_p|/m)^{k-|I|+1}\).

Let \(\Omega_0\) be the complement of the eligible boxes.
Its filtered contribution is exactly
\(\int_{\Omega_0} f(y)^m F_k(f(y)) w(y)\,dy\).
Partition at an index maximizing \(|q_i|p_i(y)\).
On that part, \(|f|\le s|q_i|p_i\) and
\(|f-z_i|\le K|z_i|\). With \(D'_i=\prod_{j\notin I}d_{ij}\)
it follows that \(|F_k(f)|\le (Km)^{s-1} D'_i d_{ip}^{\,k}\).
If \(i\) is eligible, \(p_i\le e^{-R^2}\) on its part and
(1) bounds the unfiltered integral of the modulus by
\(K e^{-c_3 m}|z_i|^m C_i m^{-\rho_i}\), absorbing \(s^m\)
and the small ratio loss \(|q_i/z_i|^m\).
If \(i\) is noneligible the same estimate with \(s^m\)
in place of \(e^{-c_3m}\), after multiplication by \(V'_i\),
is absorbed by (4).

For every \(i\) the relative factor comparison is
\[
      \frac{D'_i d_{ip}^{\,k}}{V'_i}
                  (|z_p|/m)^{|I|-1-k}=O(1).                \tag{5}
\]
Inside \(I\) the \(|I|-1\) inside factors and \(d_{ip}\) are
comparable to \(|z_p|/m\). Outside \(I\), \(D'_i/V'_i\)
equals \((|z_i|/m)/\prod_{j\in I}d_{ij}\) and each \(d_{ij}\)
there is asymptotic to \(d_{ip}\). Thus the left side is bounded
by a constant times \((|z_i|/m)/d_{ip}\) times
\(((|z_p|/m)/d_{ip})^{|I|-1-k}\).
Both base ratios are bounded by the triangle inequality and
relative separation, even when \(z_p\to0\). The filtered
contribution from \(\Omega_0\) divided by \(Z_k\) tends to zero.

For the remaining eligible terms, the exact shift formula is
\[
 F_k(E)(z_i^t A_i(t))\big|_{t=m}
   = z_i^m\sum_{h=0}^{\deg F_k}
       \frac{z_i^h F_k^{(h)}(z_i)}{h!}\Delta^h A_i(m),
                   \qquad A_i(t)=C_i t^{-\rho_i}B_i(t).
\]
For \(i\notin I\), a nonzero derivative summand must remove
its factor \(z-z_i\). Charge each derivative with a factor
\(|z_i|/m\) using (3), and bound untouched factors by their
padded distances. This gives (5) with the additional information
that \((|z_i|/m)/d_{ip}\to0\), so the normalized contribution
tends to zero. For \(i\in I\), both \(|z_i/z_p|^m\) and
\(V_i'/V_p'\) are bounded above and below. Maximality at \(p\)
bounds \(C_i m^{-\rho_i}/(C_p m^{-\rho_p})\). By the prefactor
hypothesis \(r_i>r_p\) would make that ratio divergent, and
\(r_i<r_p\) makes it tend to zero. For the surviving equal
counts pass to limits \(C_i/C_p\to\omega_i\ge0\), with
\(\omega_p=1\). Further, \(P(z_i)/P(z_p)\to1\) and any
derivative of \(P\), scaled by the corresponding power of
\(z_p/m\), is negligible compared with \(P(z_p)\), by outside
separation. Formula (3) and vanishing of the window give
\[
 \sum_{i\in I_0}\omega_i e^{v_i}
       \sum_{h=0}^k\binom{k}{h}v_i^{\,k-h}
                         (-1)^h(\rho_p)^{\overline h}=0
                       \qquad(0\le k<|I|),
\]
where \(I_0\) consists of the surviving eligible labels and
the common positive amplitude has been removed. As in the
[mixed-coordinate filter](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md#unbounded-starting-order),
the inner polynomial is monic of degree \(k\), so grouping
equal \(v_i\) gives a Vandermonde system. The group containing
\(p\) has nonzero weight \(e^{v_p}\sum_{v_i=v_p}\omega_i\).
This is a contradiction.

## Bounded starting orders

Pass instead to a fixed finite window. For a mixed count vector
\(k\) with \(\sum_i k_i=t\), the bound (1) and the variance
identity give
\[
 \prod_i p_i(y)^{k_i}\le\exp\!\left(-t\|y-\bar u\|^2
       -t^{-1}\sum_{i<j} k_i k_j\ell_{ij}^2\right),
                   \qquad \bar u=\sum_i k_i u_i/t .
\]
If \(\ell_*\) is the largest distance between participating
labels, the Gaussian integral with weight \(w\) is bounded
by \(K_t\prod_\nu(1+\bar u_\nu)\), and the latter product is
at most \(e^{r\ell_*}\prod_i C_i^{k_i/t}\). This follows by
bounding the coordinatewise log ranges of \(1+u_{\nu i}\)
by \(\ell_*\) and using \(1+u_{\nu i}\le2u_{\nu i}\) in
positive positions. The variance contributes a loss at least
\(\ell_*^2/t\). Hence all mixed terms in (1) contribute
\(o(\max_i C_i |q_i|^t)\) uniformly over the coefficients.
For pure terms direct Gamma integration and
[Stirling's formula](https://dlmf.nist.gov/5.11) give
\[
           \int p_i^t w\,dy
              =C_i t^{-\rho_i}\big((\pi/2)^{r_i/2}+o(1)\big)
\]
as the minimum positive entry tends to infinity.

Pass to limits of pairwise ratios of \(|q_i|\), \(C_i\), and
\(C_i|q_i|^t\) for every test order, allowing zero or infinity.
Group labels where both \(|q_i|\) and \(C_i\) are mutually
comparable. Among groups with mutually comparable \(|q_i|\),
discard those with \(C_i\) asymptotically smaller than another
group. They never attain a leading score \(C_i|q_i|^t\).
The remaining \(g\) groups have strictly ordered frequency
magnitude scales. If a group leads at \(t\), then at any
larger order every strictly smaller frequency scale has
score negligible compared with that group. Thus there are
at most \(g-1\) orders at which distinct groups tie for the
leading scale: each tie consumes a transition from its smallest
to its largest leading frequency scale. The strictly leading
orders of each group form an interval, by the affine dependence
of the log score ratios on \(t\).

Some group is strictly leading on at least its own cardinality
in consecutive orders. Otherwise the number of strict orders
would be at most the sum of group cardinalities minus \(g\),
which is at most \(s'-g\), whereas the window has at least
\(s'-(g-1)\) strict orders. Comparable \(C_i\) imply equal active
counts by hypothesis. At these consecutive orders normalize by
one reference \(C_p q_p^t t^{-\rho_p}\) in the leading group.
Pass to the limiting frequency ratios and the positive weight
ratios. The consecutive weighted power-sum equations violate
the Vandermonde determinant after grouping equal frequencies.
All frequencies in this group are nonzero. This contradicts
the assumed failure.

Finally, the finitely generated field of the coefficients over
\(\mathbb Q\) embeds in \(\mathbb C\), preserving vanishing of
the integer-coefficient moment forms.

## A mixed-prefactor obstruction in the leading model

If a relative cluster retains unequal \(\rho_i\), the first
three polynomial factors in its limiting filter are
\(1,w-\rho,(w-\rho)^2+\rho\). These equations can all vanish
with positive unphased weights. Take two modes \(w=\pm i\theta\)
with parameter \(\rho\) and weight one each. Take a real third
mode \(w=\beta\) with parameter \(\sigma=\rho-10\), where
\(\theta\in(\pi/2,\pi)\) satisfies \(\theta^2/\cos^2\theta=10\),
and
\[
            \beta=\sigma-\rho-\theta\tan\theta,\qquad
            \omega_3=-2\cos\theta\,e^{-\beta}>0.
\]
Such a \(\theta\) exists by continuity since the ratio tends
to infinity at the left endpoint and to \(\pi^2<10\) at the
right endpoint. The third weight cancels the zeroth sum.
The ratio of the first-order sum of the first two modes to
their zeroth sum is \(-\rho-\theta\tan\theta=\beta-\sigma\).
For the second-order sum the corresponding ratio is
\((\beta-\sigma)^2+\rho-\theta^2/\cos^2\theta
=(\beta-\sigma)^2+\sigma\).
Thus the leading three-mode system can collide. This concerns
the limiting filter equations; realizing such a collision in
ordinary factorial moments is a separate question.

These limiting scales occur on homogeneous integer supports.
In \(r=22\) variables let the first two columns each put \(d/2\)
on a pair of coordinates, with the two pairs disjoint. Let the
third column put \(A\) on \(21\) coordinates and \(d-21A\) on
the last, whose position is outside both pairs. Take \(A,m\to
\infty\) and even \(d\sim K m^{20}A^{21}\), with \(K>0\).
The powers are \(\rho=21\) for the first two columns and
\(\sigma=11\) for the third. The ratio of the third pure
saddle prefactor to either of the others tends to
\(2^{21}(\pi/2)^{10}/\sqrt K\), which can equal \(\omega_3\).
Their square-root distances are of order \(\sqrt d\).
Taking the three raw \(q_i\) to be
\(e^{i\theta/m},e^{-i\theta/m},e^{\beta/m}\) makes the
saddle shifts exponentially small on the relative
frequency scale. This produces the stated limit from
ordinary supports and coefficients, with exact zero
moments still unestablished.
