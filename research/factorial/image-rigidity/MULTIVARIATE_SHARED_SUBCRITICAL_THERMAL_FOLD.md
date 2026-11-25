# The subcritical thermal fold of a varying fully shared profile

Use the phase, parameters, and finite real integral of the
[moving quartic window](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#uniform-centers-and-finite-integral-reduction).
Consider profiles on \(0<W_-\le W\le W_+<\infty\) with
\(\delta=\max\rho-1\to0^-\) and \(\sqrt n\,|\delta|\to\infty\).
Then \(a_2<0\) and \(|a_2|\asymp|\delta|\), since
\(q_*-q_m=O(|\delta|)\) and \(\rho'(q_m)=0\).
Put
\[
 r=\sqrt{a_2/(2g)},\qquad c=-4g r^4,\qquad M=nc,\qquad
 b=\frac{r}{c}(w_*-w_n).
\]
Here \(r\to0\), \(c\asymp r^4\), and \(M\to\infty\).
The upper fold of the limiting thermal phase has field
\(b_c=2i/(3\sqrt3)\). Conjugation treats the lower fold.

**Theorem.** Along these profiles, a consecutive common root
cannot have both scaled fields \(b,b'\to b_c\), where
\(b'=r(w_*-w_{n+1})/c\) has the same imaginary part as \(b\).
The total-degree difference is arbitrary.

## A fixed cubic chart for the thermal contour

Put \(q=q_*+rz\) in the exact local phase. After removing the
constant and dividing by \(c\), it becomes
\[
 \Psi(z,b)=-z^4/4-z^2/2+bz+r z^5 E_r(z).
\]
The analytic \(E_r\) is uniformly bounded on
\(|rz|\le r_0\), with uniform derivatives on compact \(z\)-sets.
On each fixed bounded \(z\)-set the phase and its derivatives
thus converge to the displayed quartic model. The integral
over the fixed physical interval from the moving quartic
proof can be shifted to \(\Im z=1/\sqrt3\); the end
connectors have a strict physical phase loss. On that line
at \(b=b_c\) the real part of the model phase differs from
its value at \(z_c=i/\sqrt3\) by exactly \(-(\Re z)^4/4\).
At bounded heights and large \(|\Re z|\), the term
\(r z^5 E_r(z)\) is at most a small multiple of \(|z|^4\)
by the choice of \(r_0\). Thus for every sufficiently small
fixed central interval of the shifted line, the complement
has a strict scaled loss and the far tails have a quartic
majorant. These bounds persist for complex \(b\) in a small
disk about \(b_c\) and small \(r\).

The model is cubic to leading order at this fold. At
\(b=b_c\), writing \(h=z-z_c\), the stationary increment is
\[
 -ih^3/\sqrt3-h^4/4=(-i/\sqrt3)h^3(1-i\sqrt3 h/4).
\]
The canonical coordinate there is
\(u=i3^{1/6}h(1-i\sqrt3 h/4)^{1/3}\).
The endpoints at \(h=\pm s\), for any fixed small \(s>0\),
lie in the open cubic decay sectors of arguments
\(\pm\pi/3\). In fact
\(u=i3^{1/6}h+3^{1/6}\sqrt3\,h^2/12+O(h^3)\),
so both endpoints are strictly on the decay side of the
imaginary \(u\)-axis. The analytic cubic normal form
\[
                    C_r(b)+u^3/3-\lambda_r(b)u
\]
exists for the actual phase with a common chart and
parameter disk. Its coordinate and parameter converge
with derivatives to the model chart. This is the usual
two-stationary-point cubic normal form, since
\(\Psi_z=\Psi_{zz}=0\), \(\Psi_{zzz}\ne0\) and
\(\partial_b\Psi_z=1\) at the model fold. Its parameter
derivative obeys \(\lambda_0'(b_c)=i3^{-1/6}\), by matching
the linear term in \(h\). These charts have the two same
decay sectors as the model, so the local class is the
elementary Airy contour. The omitted integrals have a
strict scaled loss and can be carried as an analytic
remainder with the exact local integral.

## Uniform location and differentiation of the Airy roots

For each fixed profile interpolate the exact local integral
by replacing \(n\) with \(j=M/c\) in the amplitude of the
moving quartic proof. This gives an analytic function of
complex \(M\) in a fixed narrow sector about the positive
ray. Sectorial Stirling there keeps the amplitude uniformly
analytic and the center value bounded away from zero. At
integer \(j\) the exact moment differs from this integral by
the cotangent and excluded-summand error in that proof.

Remove the common exponential and saddle prefactor.
Put \(\xi=M^{2/3}\lambda_r(b)\), with the branch at positive \(M\).
The integral is
\[
 I(M,b)=A(M,b)\mathrm{Ai}(\xi)
        +M^{-1/3}B(M,b)\mathrm{Ai}'(\xi)+E(M,b),
 \tag{1}
\]
where \(A\) is uniformly bounded above and away from zero,
and \(B\) is uniformly bounded. Analytic division of the
amplitude by \(u^2-\lambda_r(b)\) and one integration by
parts give a remainder bounded by
\[
 |E|\le C |M|^{-1}
 \bigl(|\mathrm{Ai}(\xi)|+
       (1+|\xi|)^{-1/2}|\mathrm{Ai}'(\xi)|\bigr).
 \tag{2}
\]
The absolute-integral bound in the
[joint Airy contour proof](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#a-uniform-joint-norm-for-the-actual-cubic-class)
gives this joint norm after \(v=M^{1/3}u\) and removal of
the Jacobian. The factor \(|M|^{-1}\) is from integration
by parts. The substitution and completion of the contour
keep the ends in the two open decay sectors for a
sufficiently small fixed sector of \(M\).
The phase gap outside the chart persists in this sector:
on the horizontal tails the imaginary part of the model
is \(O(1+|\Re z|^3)\), and the perturbation is bounded
by a small multiple of the quartic loss. On bounded
contours persistence follows from the strict real-order
gap at \(b_c\). The exterior pieces and the chart
connectors are absorbed in (2) by the joint lower bound
on a smaller field disk. The lattice remainder after
normalization at integer orders is
\(\exp(-\eta n+O(M))\) times a power of \(n\), since
\(M=o(n)\), \(r\gg n^{-1/4}\), and the removed exponential
is measured relative to the quartic center.

At any actual root approaching the fold, (1) and (2)
give \(|\mathrm{Ai}(\xi)|\le C M^{-1/3}|\mathrm{Ai}'(\xi)|\).
The
[Airy continuation argument](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-elementary-airy-class-and-its-transverse-band)
puts \(\xi\) within \(O(M^{-1/3})\) of a negative real
Airy zero \(a_k\). The product guard holds because the
square of the Airy coefficient times \(|\xi|\) is
\(O(|\lambda_r(b)|)\to0\). In particular
\[
 |\Im\lambda_r(b)|=O(M^{-1}),\qquad
 |a_k|M^{-2/3}\to0 .
 \tag{3}
\]
At the next order \(M\) changes to \(M+c\).
Since \(\Im\lambda_r'\) stays strictly positive near \(b_c\)
and \(b'-b\) is real, (3) gives \(|b'-b|\le C/M\).
Thus the two \(\xi\) values differ by \(O(M^{-1/3})\).
They have the same label \(k\). Indeed
\(y(s)=\mathrm{Ai}(a_k+s)/\mathrm{Ai}'(a_k)\) solves
\(y''=(a_k+s)y,\ y(0)=0,\ y'(0)=1\). Its integral
equation gives \(y(s)=s[1+O((1+|a_k|)|s|^2)]\) on
disks of radius at most a small multiple of
\((1+|a_k|)^{-1/2}\). Distinct Airy zeros are separated
by at least that radius, which is much larger than
\(M^{-1/3}\) in (3).

At a given positive \(M_0\), fix the selected \(a_k\).
On \(|M-M_0|\le\varepsilon M_0\) make the analytic
substitution \(b=\lambda_r^{-1}(M^{-2/3}\xi)\).
This stays in the fold chart for
\(|\xi-a_k|\le R M_0^{-1/3}\), with any fixed \(R\),
by (3). The normalized Airy solution on this disk
satisfies \(y(s)=s(1+o(1))\), \(y'(s)=1+o(1)\).
Dividing (1) by \(A(M,b)\mathrm{Ai}'(a_k)\)
adds a term \(O(M_0^{-1/3})\) to \(y\).
For large fixed \(R\), Rouché gives a unique simple zero
\(\xi_k(M)=a_k+\epsilon_k(M)\), analytic in \(M\) on
the complex disk. Its \(O(M_0^{-1/3})\) displacement
and Cauchy's estimate give
\[
       |\epsilon'_k(M)|\le C M_0^{-4/3}
       \qquad (|M-M_0|\le\varepsilon M_0/2).
 \tag{4}
\]
The derivative bound uses the exact integral including
its exterior scaled pieces. Only the lattice error is
omitted. Cauchy's estimate on the same \(\xi\)-disk also
bounds the derivative of the normalized difference from
\(y\). More precisely, use a larger fixed \(R\) in (1).
The term \(M^{-1/3}(B/A)y'\) has derivative \(o(1)\)
because \((B/A)_\xi=O(M_0^{-2/3})\) and
\(y''=O((1+|a_k|)M_0^{-1/3})\). The error (2) has
derivative \(O(M_0^{-2/3})\) by Cauchy.
Thus the normalized zero slope is bounded away from
zero and the actual root at each integer order is
within \(\exp(-\eta' n)\) of the interpolating zero
in the \(\xi\)-coordinate, for some \(\eta'>0\).
Here the division by \(\mathrm{Ai}'(a_k)\) costs at
most \(\exp(O(M))\) by the lower joint bound, or by
the usual Airy derivative estimates.

## Consecutive imaginary motion

In the physical thermal field, the analytic zero branch is
\[
 b_k(M)=\lambda_r^{-1}(M^{-2/3}(a_k+\epsilon_k(M))).
\]
Equations (3) and (4) show that for real \(M\) between
\(M_0=nc\) and \(M_0+c\) its derivative is
\[
 b_k'(M)=(\lambda_r^{-1})'
    \left(-\tfrac23 a_k M^{-5/3}+O(M^{-2})\right).
\]
The inverse derivative tends to \(-i3^{1/6}\).
Since the negative Airy zeros are bounded away from zero,
the imaginary derivative has strict constant sign and
modulus at least \(c_1 |a_k| M_0^{-5/3}\)
throughout this interval. Therefore the imaginary parts
of the two interpolating roots differ by at least
\(c_1 c |a_k| M_0^{-5/3}\).
This is bounded below by a fixed positive multiple of
\(n^{-5/3}\), because \(c=M_0/n\) and \(M_0\le n\)
eventually. The exponential lattice errors cannot remove
the difference. At a common root the exact drift is
\(b'-b=(r/c)\Delta\log(1+1/n)\), which is real. The
different imaginary parts give the contradiction.
