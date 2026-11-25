# Exclusion at double zeros of the moving private cubic model

Use the ordinary one-private-linear profiles and normalized
interpolation in the [moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Write \(\epsilon=n^{-1/3}\),
\(T=\ell^2(1-\sigma)/\epsilon\), and \(H=\ell\xi_n/\epsilon^2\).
The exact order drift is real in \(\xi\).

**Theorem.** Along bounded-curvature profiles in the moving
cubic window, both scaled fields of a consecutive common root
cannot approach the same double zero of \(F(T_0,\cdot)\) at any
real limiting temperature \(T_0\). In particular, when
\(\Delta=o(n^{1/3})\) and both fields and temperatures are bounded,
a limiting common root must be simple. With \(\Delta\ne0\) in
that range there are no common roots at sufficiently large
orders. With \(\Delta=0\) a limiting common root requires
\(2H_0-T_0 H_T=0\) at a simple model zero.

## Coefficient structure at a double limit

Define the real coefficients \(c_j\) by the Taylor expansion
of the phase derivative at \(\xi=0\):
\[
 \Phi_q(q,0)=-(1-\sigma)q-c_2q^2-c_3q^3-c_4q^4+O(q^5).
\]
In particular \(c_3=(1-\sum_i d_i^4/a_i^3)/3\).
Write \({\cal J}(T,H,\epsilon)\) for the normalized
interpolation, treating \(T,H,\epsilon\) as independent
for the fixed profile. On compact sets of \(T,H\),
\[
 {\cal J}
 = F+\epsilon\{a F_H-b F_{HHHH}\}+\epsilon^2 L_2+O(\epsilon^3),
 \quad
 a=\ell A_0'(0),\quad b=c_3\ell^4/4 .                     \tag{1}
\]
The function \(L_2\) is a linear combination of derivatives
of \(F\) and the scalar \(1/(2\ell)\). Indeed substitution
\(q=\ell\epsilon y\) gives the order-\(\epsilon\) phase
term \(-b\epsilon y^4\) and order-\(\epsilon^2\) term
\(-c_4\ell^5\epsilon^2 y^5/5\). The order-\(\epsilon^2\)
amplitude is \(\ell^2 A_0''(0)y^2/2\). Their products give
the other polynomial weights at this order, and
the lattice term \(B_j/(\ell\epsilon)\) contributes
\(\epsilon^2/(2\ell)\). The order error from \(A_j-A_0\)
is \(O(j^{-1})=O(\epsilon^3)\).

These expansions have uniform derivatives in \(T,H\) on
smaller compact sets, and retain two bounded derivatives
in real \(\epsilon\). For the integral part, the real
segment and its differentiated Taylor integrands are
dominated by a polynomial times
\(\exp(-c y^3+C y^2+C|y|)\), by the same fixed endpoint
cutoff and sectorial Stirling estimates as in (5) of the
moving model proof. The analytic order derivatives of
\(A_j-A_0\) follow by Cauchy from its \(O(j^{-1})\) bound.
The normalized lattice correction is the
explicit \(B_j/(\ell\epsilon)\) in that proof. Its remainder
after \(\epsilon^2/(2\ell)\) is \(O(\epsilon^4)\) on
complex \(\epsilon\) disks of radius proportional to real
\(\epsilon\), by (5) with independently substituted
linear and quadratic phase terms. Cauchy estimates also
control its first two derivatives. Thus Taylor expansion
and differentiation may be used jointly in (1).

There are two types of ordinary profiles near \(\sigma=1\).
If the shared row \((1,2)\) is absent, all shared
\(a_i\ge2\) once \(\sigma\) is sufficiently near \(1\).
Indeed a distinct positive row with \(a_i=1\) and
\(d_i\ne1\) contributes at least \(4\) to \(\sigma\).
Thus \(\sum_i d_i^4/a_i^3\le\sigma^2/2\) and \(c_3\)
stays strictly positive.
If \((1,2)\) is present, it is unique and
\(u=\sigma-1\ge0\) is the extra shared curvature.
The remaining sums satisfy
\(\sum d_i^4/a_i^3\le u^2\),
\(\sum |d_i|^3/a_i^2\le u^2\), and
\(\sum |d_i|/a_i\le u\) on the extra rows.
Consequently \(c_3=O(u^2)\), \(c_2\to1\), \(T\le0\),
\(a\to1\), and the coefficients of \(L_2\) tend to those
from \((1,2)\) and the private linear row alone:
\[
                         L_2\longrightarrow
                F_{HH}/2-F^{(5)}/10+1/2 .                 \tag{2}
\]
Here \(c_4=1/2\) and \(A_0(q)=\sqrt{(1+q)/(1-q)}\) for
those two rows, giving (2). Along the bounded window
\(u=O(\epsilon)\), so the phase part \(b\) of the
linear correction in (1) is \(O(\epsilon^2)\).

## Real degree drift and the quadratic discriminant

Fix a double zero \((T_0,H_0)\), let \(t=T-T_0\) and
\(z=H-H_0\), and orient by conjugation so \(\Im H_0>0\).
The identity \(F_{HH}+T F_H-HF=1\) gives
\(F_{HH}=1\) and \(F_T=-1/2\) at the double zero.
A small fixed circle about \(H_0\) encloses the two roots
of the local interpolation, counted with multiplicity.
Let \({\cal W}(z)=z^2+A z+B\) be their monic polynomial
and \({\cal D}=A^2-4B\). The coefficients and their
parameter derivatives are uniformly bounded by
integrating \({\cal J}_H/{\cal J}\) and its derivatives
on that circle. The interpolating factor
\({\cal J}/{\cal W}\) stays bounded away from zero on
a smaller disk, since its limit at the double zero
is \(1/2\).

At the two orders the coefficients of \({\cal W}\)
change by \(O(\epsilon^3)\), since \(T_N-T_n=O(1/n)\)
and \(\epsilon_N-\epsilon_n=O(\epsilon/n)\).
Write
\[
 H_N-H_n=v+O(\epsilon^3),\qquad
 v=\ell N^{2/3}\Delta\log(N/n)\in\mathbb R .
\]
If \(\Delta\ne0\), then \(|v|\gtrsim\epsilon\) and, when
both fields approach \(H_0\), \(v\to0\).
At an actual common root both corresponding polynomial
values are exponentially small. Subtract them using
the lower-order polynomial at both arguments. This gives
\[
 2z+A=-(H_N-H_n)+O(\epsilon^3/|v|),
 \qquad
 {\cal D}=v^2+O(\epsilon^3).                               \tag{3}
\]
In fact \((2z+A)^2-{\cal D}=4{\cal W}(z)\). Thus (3)
uses the full linked change of orders even if \(v\)
tends to zero arbitrarily slowly.

For the limiting model, the two local branches of the
double zero have expansions
\[
 z=\pm \sqrt t-T_0t/3
       \pm (H_0/12+T_0^2/36)t^{3/2}+O(t^2).
\]
These follow by substituting into the Taylor expansion
of \(F\) using \(F_T=-F_{HH}/2\) and the differential
equation. Equivalently they are the
[double-zero branches](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#isolated-exceptional-temperatures).
Squaring their difference gives the analytic model
discriminant
\[
 {\cal D}_0(t)
 =4t+ (2H_0/3+2T_0^2/9)t^2+O(t^3).
 \tag{4}
\]
To find the first perturbation at \(t=0\), the unit
in \({\cal J}/{\cal W}\) has value \(1/2\) at the
base double point. Thus the perturbation of
\({\cal D}\) is \(-8\) times the perturbation of
\({\cal J}(T_0,H_0)\) to first order. In (1),
\(F_H=0\) and \(F_{HHHH}=H_0+T_0^2\) at the double
zero. It follows that, for profiles without \((1,2)\),
\[
 \Im{\cal D}=(\Im H_0)\left(\tfrac23 t^2
                     +2 c_3\ell^4\epsilon\right)
                  +O(|t|^3+\epsilon|t|+\epsilon^2).
 \tag{5}
\]
The leading terms have the same sign. Since \(c_3\)
and \(\ell\) stay bounded away from zero, (5)
contradicts (3).

In profiles containing \((1,2)\), translating the field
by \(a\epsilon\) absorbs \(\epsilon a F_H\) in (1) and
preserves the discriminant. The remaining first-order
phase term is \(O(\epsilon^3)\) on the bounded window.
Thus in these profiles
\[
 {\cal D}(t,\epsilon)
 ={\cal D}_0(t)
      +\epsilon^2(-8/5-4T_0^3/5-8T_0 H_0/5)
      +o(\epsilon^2)
 \tag{6}
\]
along \(t\to0\). To compute the coefficient, the translated
interpolation differs from \(F\) to second order by
\(L_2-a^2 F_{HH}/2\); the limit (2) makes its value at
the double zero \(1/2-F^{(5)}/10\).
The equation for \(F\) gives
\(F^{(5)}=3-T_0^3-2T_0 H_0\) there.
Multiplication by \(-8\) yields (6).
Because \(T_0\le0\), when \(T_0<0\) the two leading
imaginary terms in (4) and (6) have the same sign,
again contradicting (3). When \(T_0=0\), one has
\(t=T\le0\); the real part of (6) reads
\(4t+O(t^2)-8\epsilon^2/5+o(\epsilon^2)\),
which is strictly negative on the scale \(|t|+\epsilon^2\).
This also contradicts (3).

## The equal-degree increment

If \(\Delta=0\), the physical field \(\xi_j\) is the same
at the two orders. Thus \(T_j\), \(H_j\), and
\(\epsilon_j\) change by leading relative multiples
\(1/(3n)\), \(2/(3n)\), and \(-1/(3n)\), respectively.
Subtracting the two interpolated values, multiplying
by \(3n\), and using the bounded second derivatives
from (1) gives
\[
 T{\cal J}_T+2H{\cal J}_H-\epsilon{\cal J}_\epsilon
                                      =O(\epsilon^3).     \tag{7}
\]
The exponentially small integer-order errors are
absorbed here. At a double limit \(F_H=0\),
\(F_T=-1/2\), so (7) forces \(T_0=0\).

Now \(F_{HH}=1\) and \(F_{HT}=0\) at \((0,H_0)\).
By (1) and (7), \(z=O(|T|+\epsilon)\), since
\(H_0\ne0\). Substitution in
\({\cal J}=O(e^{-\eta n})\) first gives
\(|T|=O(\epsilon)\).
For a profile without \((1,2)\), this reduces the zero
equation to
\(-T/2-\epsilon b H_0=O(\epsilon^2)\).
Its imaginary part is impossible because \(b\) stays
strictly positive and \(\Im H_0>0\).

For a profile containing \((1,2)\), the order-\(\epsilon\)
phase correction is \(O(\epsilon^3)\), and the
linear amplitude term is \(O(\epsilon(|z|+|T|))\).
The zero equation now improves the bound to
\(|T|=O(\epsilon^2)\). Equation (7) then gives
\(z=-a\epsilon+o(\epsilon)\), since its leading term
is \(2H_0(z+a\epsilon)\).
At \(T_0=0\), (2) has value \(7/10\) at the double zero.
Thus the original zero equation becomes
\[
                  0=-T/2+\epsilon^2/5+o(\epsilon^2),
\]
also impossible since \(T\le0\).

If \(\Delta=o(n^{1/3})\) and both fields are bounded, the
two fields necessarily tend to the same zero along a
subsequence. The theorem leaves only a simple limit.
There the analytic zero branch of the interpolation
changes the physical field by \(O(n^{-5/3})\), as proved
in the moving cubic model. This excludes nonzero
\(\Delta\), whose exact real drift has modulus at least
\(\log(1+1/n)\). For \(\Delta=0\), the limit of (7) is
\(T_0 F_T+2H_0 F_H=0\). At a simple zero this is
\(2H_0-T_0 H_T=0\), as asserted.
