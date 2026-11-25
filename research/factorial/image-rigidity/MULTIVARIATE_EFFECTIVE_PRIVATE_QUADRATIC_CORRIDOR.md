# An effective full private quadratic corridor

Use the one-private-linear profile, common physical field
convention, and \(\mu=1-\sigma\) of the
[bounded quadratic exclusion](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_WINDOW.md).
Suppose \(0<\mu\le1\) and \(\xi_j\ne0\).
Put
\[
    r_j=\frac{|\xi_j|}{\mu},\qquad
    M_j=j\mu r_j^2,\qquad H_j=\sqrt{j/\mu}\,\xi_j ,
\]
so \(M_j=|H_j|^2\).

**Theorem.** There are computable positive rational \(p_*,K\)
such that no consecutive common root occurs when
\[
             |\xi_j|/\mu^2\le p_*\quad(j=n,N),
                     \qquad n\mu^3\ge K.                   \tag{1}
\]
There is no restriction on the integer degree difference.
This includes both bounded and escaping quadratic fields.
It leaves the thermally scaled fields with \(|\xi|/\mu^2\)
bounded away from zero outside this corridor.

## A computable escaping threshold

We first require \(M_j\ge M_0\) at both orders for a computable
\(M_0\), in place of \(n\mu^3\ge K\).
Use separate small bounds \(r_j\le r_0\) and \(r_j/\mu\le p_0\)
during the proof. Both follow from (1) on taking
\(p_*\le\min(r_0,p_0)\), since \(0<\mu\le1\).
Use the fixed \(L,q_0\), phase and Gamma amplitude bounds of the
[bounded quadratic proof](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_WINDOW.md#a-real-order-interpolation-with-a-fourth-order-lattice-error).
In particular \(1/4\le c_2\le1\), \(|E|\le8\) on \(|q|\le1/4\),
and \(\Re\Phi\le-\gamma<0\) throughout the real tail \(q\ge q_0\)
for small physical fields.

Use the effective finite cotangent representation (2)–(3) of the
[private cubic field proof](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#effective-lattice-representation-on-a-fixed-small-rectangle).
It holds with the same choice rules on the quadratic corridor.
On the real segment \([0,L]\) the quadratic term is nonpositive,
and so are the cubic and higher terms combined. Thus the upper
bound on the real phase maximum used there can be taken as
\(|\xi|L\) without requiring \(\mu\) small.
The right and horizontal wall checks only use the uniform bound
on \(\Phi_{qq}\), the real loss \(\gamma\), and a sufficiently
small physical field, all available here.
Consequently the normalized exact moment is
\(\int_0^L A_j(q)e^{j\Phi(q,\xi)}\,dq+B_j(\xi)\)
up to \(C e^{-\gamma_1j}\), with bounded amplitude and its
order derivative \(O(j^{-2})\) and with the same effective
bounds \(|B_j|,|\partial_\xi B_j|\le C/j\),
\(|\partial_j B_j|\le C/j^2\).

At an individual order use its own scales
\(r=r_j,\ c=\mu r^2,\ M=jc,\ b=\xi/|\xi|\).
With \(p=c_2 r/\mu\), the rescaled phase is
\[
     \Psi_b(z)=bz-z^2/2-p z^3/3+(rp/c_2)z^4 E(rz).
\]
The profile hypotheses make \(p\le p_0\).
After conjugation \(\Im b>0\), since real \(\xi\) in the field
neighborhood corresponds to positive \(t\).
The pure quadratic phase has saddle \(z_0=b\).

If \(\Re b\le-1/4\), the real integration interval gives a
strict linearly descending endpoint contour. The last two phase
terms combined contribute nonpositively on the positive real
axis, by \(8L/c_2<1/12\). Thus the endpoint term is
\(-r/(M b)(1+O(1/M+c/r))\), with the lattice correction included.
For \(-1/4\le\Re b\le1/4\), move to the horizontal line at
height \(1/2\), starting along the vertical leg from zero.
Here \(\Im b\ge\sqrt{15}/4\).
The pure vertical phase \(-v\Im b+v^2/2\) has strict linear
loss for \(0<v\le1/2\), and the pure horizontal real phase
is bounded above by
\(1/32-\sqrt{15}/8+1/8<0\).
Both bounds persist for small \(p_0,r_0\), using the displayed
phase and its derivatives on a bounded scaled set.
The tail bound below supplies the same conclusion on the rest
of the horizontal line. This again gives only the nonzero endpoint
term, up to a relative error tending effectively to zero.

For \(1/4\le\Re b\le1\), Rouché around the pure saddle
continues a simple saddle \(z_s=x+iy\) with
\(z_s=b+O(p)\), and \(\Re(-\Psi_b''(z_s))\) bounded below.
For example use radius \(1/40\) around \(b\), where the
unperturbed phase derivative has modulus \(1/40\).
The derivative perturbation is bounded by \(C p\) there
using the fixed analytic bounds.
For positive real \(b\) this root is real and moves initially
into the upper half-plane as \(\Im b\) increases, by
\(\partial_b z_s=-1/\Psi_b''(z_s)\).
It cannot reach the real axis while \(\Im b>0\).
Deform along \(0\to iy\) and the rightward horizontal line
through \(z_s\). The right connector has fixed physical
exponential loss by a sufficiently small choice of \(r_0\).
The pure horizontal real loss from its saddle is
\((s-\Re b)^2/2\), which gives strict effective loss outside
a fixed saddle interval on \(0\le s\le16\).
For \(s\ge16\) and \(rs\le L\) one has
\[
               \Re\Psi_b(s+iy)\le-s^2/4-p s^3/8
                                                               \tag{2}
\]
on all these horizontal lines, including the endpoint-only line.
Indeed \(|y|\le2\), \(|s+iy|\le9s/8\),
\(\Re(s+iy)^2\ge63s^2/64\), and
\(\Re(s+iy)^3\ge15s^3/16\).
Use \(|b|\le2\), \(r|s+iy|\le9L/8\), and \(c_2\ge1/4\)
to bound the linear and quartic terms in (2).
Thus all Gaussian localization errors are effective.

When \(\Re b\ge7/8\) on the unit arc, the pure saddle value
\(\Re(b^2/2)\) exceeds \(1/4\).
On the vertical leg the pure real phase is nonpositive
as long as its height is at most \(\Im b\).
Moving its endpoint by \(O(p)\) and restoring the perturbation
changes that upper bound by \(O(p)\).
Consequently the active saddle contribution dominates and
there is no zero in this region.
For \(1/4\le\Re b\le7/8\) the height \(\Im b\) is uniformly
positive and the vertical leg has strict linear loss for small
\(p_0,r_0\). This follows from the negative linear term
\(-v\Im b\), the bound \(v\le y_s=\Im b+O(p)\), and the
bounded phase perturbation divided by \(v\).
In this region the endpoint and saddle terms are respectively
\[
   \frac rM e(b)(1+O(1/M+c/r)),\qquad
   \frac r{\sqrt M}s(b)e^{M D(b)}(1+O(M^{-1/2}+1/j)),
\]
where \(e(b)=-1/b\),
\(D(b)=\Psi_b(z_s)\), and
\(s(b)=A_0(r z_s)\sqrt{2\pi/(-\Psi_b''(z_s))}\).
Their uniform effective errors follow from the same
linear and Gaussian moment bounds as in the cubic field proof.
In particular at a zero
\[
 M\Re D(b)=-\tfrac12\log M+\log|e(b)/s(b)|
                        +O(M^{-1/2}+c/r).
\]
Since \(D(b)=b^2/2+O(p)\), the cosine of twice the angle
confines every such zero to within \(1/100\) of the ray
\(\arg b=\pi/4\), by an effective choice of thresholds.
This also proves the single-order confinement asserted below.

## Consecutive comparison and the mixed scales

When both orders have \(M_j\ge M_0\), their equal physical
imaginary coordinates put both scaled fields in a disk of
radius \(1/8\) around \(e^{i\pi/4}\) on the common lower-order
scale \(r,c=\mu r^2,M=nc\).
The preceding losses and amplitude bounds hold throughout
this disk by the same inequalities with
\(\Im b\) and \(\Re b\) uniformly positive.
In particular \(\Re D'(b)=\Re z_s(b)\) is bounded below.
The
[linked endpoint and saddle comparison](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#linked-endpoint-and-saddle-ratios)
now applies with \(d=rN\Delta\log(N/n)\).
Here are the scale checks for its errors.
Subtracting the two real cancellation equations controls
\(|d|\) by a computable multiple of
\(M^{-1/2}+c/r+(1+\log M)/n\), after absorbing the
Lipschitz amplitude term \(C|d|/M\).
The quantity \(d\) can therefore be kept bounded by one.
On the fixed lower-order contours the exact multiplier is
\(\exp(c\Psi_b(z)+dz)\) as before.
The endpoint and saddle errors are the same as (7) and (8)
in that proof. The only change on the long horizontal segment
is that \(|\Psi_b|\le C(s^2+p s^3)\), which is absorbed by
the loss (2) in the same integral estimates.

For nonzero \(\Delta\) one uses \(|d|\ge r\).
The potentially scale-sensitive quotients there satisfy
\[
   c/r^2=\mu\le1,\qquad c/r=\mu r,\qquad
   \frac{n^{-2}}r=\frac{\mu^2r^3}{M^2},\qquad
   \frac{M}{r^2}e^{-\gamma_1n}=n\mu e^{-\gamma_1 n}.
\]
Thus all the ratio errors divided by \(|d|\) are uniformly
small for large \(M_0\) and small \(r_0\), as is
\(|c\Re D|/|d|\). The saddle exponent \(d\Re z_s\)
cannot be canceled. For \(\Delta=0\) all ratio errors are
\(O(1/n)\), while
\(c\Re D=-(\log M)/(2n)+O(1/n)\) gives a contradiction.
This proves the escaping exclusion with effective \(r_0,p_0,M_0\).

Finally suppose (1) holds and at least one order has
\(M_j<M_0\). If the other order is above \(M_0\), the
angle confinement gives \(\Im H_{j'}\ge |H_{j'}|/2\)
after conjugation. The ratio of the two imaginary scaled
fields is \(\sqrt{j'/j}\le\sqrt2\), since the physical
imaginary coordinates agree. Thus both fields have
\(|H_j|\le3\sqrt{M_0}\). This same cap holds when both orders
are below \(M_0\).
Choose a rational \(B>3\sqrt{M_0}\) and then take \(K=K(B)\)
from the bounded quadratic theorem.
It excludes every such mixed or bounded pair under (1).
The constants are computed in this order, so the choice
of the bounded window cap introduces no circular dependence.
