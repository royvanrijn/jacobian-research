# Effective exclusion in the private cubic field corridor

Orient ordinary profiles with one private linear row as in the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Use \(\mu=1-\sigma\), \(c_2\), \(\xi_j\), and the same local
logarithm for both orders. Suppose \(3/4\le\sigma\le5/4\).
For either order with \(\xi_j\ne0\), set
\[
   r_j=\sqrt{|\xi_j|/c_2},\qquad a_j=\frac{\mu}{c_2r_j},
                    \qquad M_j=j c_2r_j^3 .
\]

**Theorem.** There are computable positive rational \(r_0,a_0,M_0\)
such that the two moments cannot share a physical root if both
orders satisfy
\[
                 r_j\le r_0,\qquad |a_j|\le a_0,\qquad M_j\ge M_0.
                                                               \tag{1}
\]
The degree difference can be any integer. This is the effective
cubic-field part of the
[endpoint thermal-tail comparison](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md).
The conditions require both physical fields to stay near the
private endpoint. They do not locate those fields for general
one-private profiles.

## Effective lattice representation on a fixed small rectangle

Use \(L=1/1024\) and the phase, Stirling amplitude, and explicit
Taylor bounds of the
[bounded cubic calculation](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md).
In particular \(7/64\le c_2\le41/32\), \(|E(q)|\le8\)
on \(|q|\le1/4\), and
\[
 \Phi(q,\xi)=\xi q-\mu q^2/2-c_2q^3/3+q^4E(q).
\]
For \(|\mu|\) and \(|\xi|\) sufficiently small, the real part of the
phase is at most \(-\gamma<0\) for \(q\in[L/2,1]\), with the
explicit loss of that calculation.
The bounds in (1) make the two physical parameters arbitrarily
small by reducing \(r_0\).

Use the finite cotangent contour of
[the moving model, equation (5)](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md#from-the-lattice-to-a-linked-analytic-integral).
For fixed sufficiently small rational \(\nu>0\) it gives, at
integer \(j\), the normalized exact moment
\[
   I_j(\xi)=\int_0^L A_j(q)e^{j\Phi(q,\xi)}\,dq
                         +B_j(\xi)+{\cal E}_j(\xi),
 \qquad |{\cal E}_j|\le C e^{-\gamma_1 j} .                 \tag{2}
\]
The constants in this representation can be chosen effectively:
on the rectangle the positive real parts of all Gamma arguments
are at least \(j/2\), or \(ja_i/2\) for shared row \(i\).
The derivative estimates from the real Binet formula also hold
with \(\Re x\) in place of \(x\) for \(\Re x>0\), by integrating
the reciprocal powers on the same vertical segment.
Thus the exact amplitude is uniformly bounded in modulus and
away from zero on the rectangle, obeys \(A_j(0)=1\),
and has real-order derivative \(O(|q|j^{-2})\).
The factor \(|q|\) follows by subtracting at \(q=0\) in each
Stirling remainder derivative and integrating along the segment
from zero to \(q\), as in the
[quadratic interpolation](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_WINDOW.md#order-derivative-and-model-residuals).

For clarity, let \(C_0\) bound the second phase derivative there,
computed by the logarithm series. Relative to a real point \(x\),
the real phase on its vertical fiber increases by at most
\(|\xi|\nu+C_0\nu^2/2\): the derivatives of the zero-field phase
on the real line are real. Choose \(\nu\) so that this increase is
less than \(\gamma/2\) on the right wall and less than \(\pi\nu/2\)
on the horizontal sides, by also reducing the field cap if needed.
Next require \(|\xi|L+|\mu|L^2/2<\pi\nu/2\), which bounds the
maximum of the real phase on \([0,L]\).
The nonconstant cotangent terms on the horizontal sides then
have a loss at least \(\pi\nu\) in the exponent, while the right
wall retains a loss at least \(\gamma/2\).
On that wall the cotangent kernel is its constant term plus
\(O((e^{2\pi j y}-1)^{-1})\) above and below the real axis.
The upper and lower nonconstant terms pair with opposite signs.
Their numerator difference is \(O(jy)e^{-\gamma j/2}\) for
small \(y\), by the derivative bound on \(f_j=A_j e^{j\Phi}\).
Thus the apparent pole at \(y=0\) causes a bounded error term.
Residues outside \(L\) obey the real lattice tail bound.
These comparisons give explicit constants for (2), including
the allowed displacement of the right wall to a mesh point
within \(1/j\) of \(L\).

The retained left-wall correction satisfies computable
\[
     |B_j|\le C/j,\qquad
     |\partial_\xi B_j|\le C/j,\qquad
     |\partial_jB_j|\le C/j^2 .                            \tag{3}
\]
Indeed its integral is
\(i\int_0^\nu [f_j(iy)-f_j(-iy)]/(e^{2\pi j y}-1)\,dy\)
with the half endpoint \(1/(2j)\) added.
On this segment the numerator is bounded by
\(Cy(1+j(|\xi|+y))e^{j y}\), by the derivative of \(f_j\)
and a small fixed cap on the field and height.
For \(\partial_\xi\), differentiating the two exponentials
gives \(C j y e^{jy}\).
For \(\partial_j\) in the numerator one has
\(C (|\xi|y+y^2+y/j^2)e^{jy}\), using the amplitude derivative
and the phase expansion at zero. Differentiating the denominator
gives an extra factor
\(2\pi y e^{2\pi j y}/(e^{2\pi j y}-1)^2\).
In each bound substitute \(u=jy\); the integrable majorants
near zero and infinity are elementary and give (3).
This derivation uses real-order differentiation only.

## Angular confinement of each root

Work first at one integer order and use its own
\(r=r_j,\ c=c_2r^3,\ M=jc,\ b=\xi_j/(c_2r^2)\), so \(|b|=1\).
After conjugation take \(\Im b>0\); real physical fields in this
small neighborhood correspond to positive real \(t\).
In the scaled variable \(z=q/r\),
\[
       \Psi_b(z)=bz-\tfrac12 a z^2-z^3/3
                      +(r/c_2)z^4E(rz),\qquad a=a_j .
\]
All compact bounds below are uniform as \(a,r\to0\).

When \(\Re b\le-1/4\), the real line is a strict endpoint
contour: for \(s\ge0\), \(rs\le L\),
\[
 \Re\Psi_b(s)\le -s/4+|a|s^2/2-s^3/4\le-15s/64
                    \qquad(|a|\le1/4).
\]
Here \(8L/c_2<1/12\), and complete the square in \(s\).
Watson's leading endpoint integral is \(-r/(M b)\)
with relative error at most \(C/M\), by the linear loss
and \(\Psi_b(z)-bz=O(z^2)\) near zero.
The additive error \(B_j\) is relatively \(O(c/r)\);
the physical tail in (2) is relatively \(O((M/r)e^{-\gamma_1j})\).
Thus no zero lies in this region for small \(r_0\) and large \(M_0\).

Otherwise the pure cubic saddle \(z_0=\sqrt b=x_0+iy_0\)
has \(x_0\ge\sqrt{3/8}\). The actual phase has a continued
simple saddle \(z_s=x+iy\) close to \(z_0\), with
\(x\) bounded below, \(0<y<2\), and \(\Re(-\Psi_b''(z_s))\)
bounded below. There are computable uniform bounds here,
for instance Rouché on \(|z-z_0|=1/40\) uses
\(|b-z^2|\ge(2|z_0|-1/40)/40\) and the uniform \(O(|a|+r)\)
perturbation of the derivative. The estimate
\(z_s-z_0=O(|a|+r)\) follows from the factored equation
\(b-z_s^2=O(|a|+r)\).
For positive real \(b\) the continued saddle is real.
Along the upper \(b\) domain it cannot cross the real axis,
since the phase derivative at real \(z\) has imaginary part
\(\Im b\); near the positive real ray \(\partial_b z_s\) is
positive real. Thus \(y>0\) for \(\Im b>0\).
The conclusions hold on compact neighborhoods of the indicated
arc, including a complex neighborhood of \(e^{i\pi/3}\).

Deform the real interval to the vertical leg \(0\to iy\),
the horizontal line through \(z_s\) to \(L/r+iy\),
and the right connector back to \(L/r\).
The connector has exponential physical loss by (2), since
its physical height is at most \(2r\) and can be kept below \(\nu\).
On the pure horizontal line through \(z_0\), the real loss from
the saddle at \(s+iy_0\) is
\[
                         (s-x_0)^2(s+2x_0)/3 .
\]
This strictly positive loss away from the saddle persists
on the actual line on \(0\le s\le16\), by the uniform
\(O(|a|+r)\) perturbations and the saddle's bounded nonzero Hessian.
For \(s\ge16\) one has \(\Re\Psi_b(s+iy)\le-s^3/8\)
through the cutoff: use \(|s+iy|\le9s/8\),
\(\Re(s+iy)^3\ge15s^3/16\), \(|b|\le2\),
\(|a|\le1/4\), and \(r|s+iy|\le9L/8\) in the quartic bound.
These estimates are uniform also in the stated neighborhood
of the equality ray.

Write \(D(b)=\Psi_b(z_s)\). Uniform quadratic localization
on this horizontal line gives its saddle integral
\[
        \frac r{\sqrt M}\,s(b)e^{M D(b)}
                           (1+O(M^{-1/2}+1/j)),\qquad
 s(b)=A_0(r z_s)\sqrt{2\pi/(-\Psi_b''(z_s))}.
                                                               \tag{4}
\]
All constants are computable using the preceding losses and the
phase derivatives on a bounded \(z\)-set. For example on a fixed
saddle interval bound the cubic Taylor remainder by \(C|s-x|^3\),
choose the interval so that this is bounded by a fixed fraction
of the quadratic loss, and integrate the error using Gaussian
moments. The complementary pieces have fixed loss or the cubic
tail. The same argument gives an absolute integral upper bound
on the scale of (4). The amplitude and its reciprocal are bounded.

If \(\Re b\ge3/4\), the model value
\(\Re(2b^{3/2}/3)\) has a positive minimum on the upper arc.
The real exponent on the vertical leg is at most
\(O(|a|+r)\), since \(\Re(biv)\le0\) there and \(0\le v\le y<2\).
Thus (4) strictly dominates the rest and there is no zero.
In the remaining region \(-1/4\le\Re b\le3/4\),
\(\Im b\) is bounded below, and the vertical leg has
strict linear loss by reducing \(a_0,r_0\).
Its length is bounded below. The endpoint integral,
including \(B_j\), has the form
\[
        \frac rM\, e(b)(1+O(1/M+c/r)),\qquad e(b)=-1/b .
                                                               \tag{5}
\]
The formula follows by expanding the phase to its linear term
at zero on the leg; the curvature error integrates to \(O(M^{-2})\)
and the amplitude equals one at zero with uniformly bounded
derivative. Both facts can be bounded by exponential moments
after taking absolute values. The term (3) has the relative
size \(c/r\).
Consequently a zero satisfies
\[
 M\Re D(b)=-\tfrac12\log M+\log|e(b)/s(b)|
                        +O(M^{-1/2}+c/r).                 \tag{6}
\]
The residual physical tail in (2) is absorbed here since
\(M/r\le Cj\).
On this arc \(|D(b)-2b^{3/2}/3|\le C(|a|+r)\); the real
part of \(2b^{3/2}/3\) has a unique zero at \(\arg b=\pi/3\).
Thus (6) confines every zero to any fixed small angular
neighborhood of \(\pi/3\), with computable thresholds.
For example one may use width \(1/100\) and bound the
minimum absolute real model value on the excluded closed arcs.
It is positive and can be bounded below by bisection or directly
by the displayed elementary cosine.

## Linked endpoint and saddle ratios

At a common root use the lower-order scale \(r,c,M=nc\)
at both orders. The equal imaginary physical fields and the
angular confinement place \(b=\xi_n/(c_2r^2)\) and
\(b'=\xi_N/(c_2r^2)\) in a fixed small disk around \(e^{i\pi/3}\).
Indeed their modulus ratio is the ratio of the sines of the
two angles. With the width \(1/100\) this disk can be taken
of radius \(1/8\). The continuation, strict vertical and
horizontal losses, and (4)–(6) hold uniformly throughout
this disk, by the same inequalities for the pure saddle, with
the parameter cuts separated from the disk.
Use \(M'=Nc\). The saddle height derivative is \(D'(b)=z_s(b)\),
so its real part is bounded below by an effective \(\kappa>0\)
throughout the disk. Also \(\log|e(b)/s(b)|\) is uniformly
Lipschitz, using the analytic implicit derivative and the
nonzero denominators in (4) and (5).

Put \(d=M'(b'-b)=rN\Delta\log(N/n)\in\mathbb R\).
Subtracting (6) at the two orders gives
\[
 \kappa |d|\le
 C\bigl((1+\log M)/n + |d|/M +M^{-1/2}+c/r\bigr).
\]
Indeed monotonicity along the real field segment controls
\(M'\Re(D(b')-D(b))\). The remaining excess exponent
\(c\Re D(b)\) is \(O((1+\log M)/n)\) by (6).
The right-hand side uses the Lipschitz bound and the
two bounded relative errors. Since \(M\le n\) and \(c/r=O(r^2)\),
we can arrange \(|d|\le1\) by shrinking \(r_0\) and enlarging \(M_0\).

Now take the lower-order vertical leg and horizontal contour
at both orders, and call the exact split integrals, with \(B_j\)
added on the left, \(E_n,R_n\) and \(E_N,R_N\).
Apart from an amplitude ratio \(1+O(n^{-2})\), the upper
integrands differ from the lower by the exact factor
\(\exp(c\Psi_b(z)+dz)\).
On the vertical leg the increment is \(O((c+|d|)|z|)\),
so using the linear loss and (3), (5) gives
\[
 \frac{E_N}{E_n}
 =1+O\!\left(\frac{c+|d|}{M}+n^{-2}
                  +\frac{c}{rn}+\frac{c|d|}{r^2n}\right).
                                                               \tag{7}
\]
In the \(B_j\) difference use \(|\xi_N-\xi_n|=|d|/(rN)\)
and divide the resulting bound by the endpoint scale \(r/M\).

Near the saddle the increment minus its central value is
\(O(c|z-z_s|^2+|d||z-z_s|)\), so the Gaussian bound and (4) give
\[
 \frac{R_N}{R_n}
   =e^{cD(b)+d z_s(b)}
        \left(1+O\!\left(\frac cM+\frac{|d|}{\sqrt M}+n^{-2}\right)\right).
                                                               \tag{8}
\]
For (7) and (8) apply \(|e^u-1|\le |u|e^{|u|}\) and
integrate the explicit linear and Gaussian majorants.
On bounded complementary horizontal pieces the fixed phase gap
gives the same error bound after subtracting the central
multiplier. On the unbounded part \(s\ge16\), the scaled
\(|\Psi_b|\) and its difference from the saddle value are
\(O(s^3)\), by the quartic bound used above. The increment
is bounded by \(C(c s^3+|d|s)\), which the original
loss \(M s^3/8\) absorbs for \(M\) large.
Integration after subtracting the multiplier retains a factor
of \(c\) or \(|d|\). This gives the stated effective
bounds even when \(|d|\) is very small.

If both moments are zero, the two ratios agree up to an error
\(O((M/r)e^{-\gamma_1 n})\), by (2) and comparability
of \(|R_n|\) with \(|E_n|\) at the lower zero.
If \(\Delta\ne0\), then \(|d|\ge r\).
Divide the errors in (7), (8), and the exponential remainder
by \(|d|\). They can be made arbitrarily small uniformly by
increasing \(M_0\) and shrinking \(r_0\).
For example
\(n^{-2}/|d|\le C r^5/M^2\),
\(c/(rn|d|)\le c_2r/n\), and the exponential
remainder after division is at most \(C n r e^{-\gamma_1n}\).
Also \(|c\Re D(b)|/|d|\le Cr^2\).
Hence the real logarithm of the modulus of (8) has leading
size at least \(\kappa|d|\) and cannot equal that of (7).
If \(\Delta=0\), then \(d=0\) and (7), (8) and the exponential
remainder have \(O(1/n)\) errors. By (6),
\[
             c\Re D(b)=-\frac{\log M}{2n}+O(1/n).
\]
Increasing \(M_0\) until \(\log M\) strictly exceeds twice
the sum of these effective error constants gives a contradiction.
This proves the claimed computable threshold.
