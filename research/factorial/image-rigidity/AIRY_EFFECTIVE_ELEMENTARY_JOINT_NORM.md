# Effective absolute bounds on an elementary Airy contour

Write
\[
       p_z(v)=v^3/3-zv,\qquad
       J(z)=|\operatorname{Ai}(z)|
                   +(1+|z|)^{-1/2}|\operatorname{Ai}'(z)|.
\]
Let \(C_H\) be the standard truncated elementary Airy chain
from \(H e^{-i\pi/3}\) to zero to \(H e^{i\pi/3}\).
The same estimates hold for rotations of the chain by a cube
root of unity, after rotating the field.

**Lemma.** There are computable positive constants \(C,\delta,H_0\)
such that for \(H\ge H_0,\ |z|\le\delta H^2\), a function \(f\)
holomorphic on a neighborhood of \(|v|\le H\) satisfies
\[
      \left|\int_{C_H} f(v)e^{p_z(v)}\,dv\right|
             \le C\sup_{|v|\le H}|f(v)|\,J(z).              \tag{1}
\]
The same bound holds with the endpoints' angles rotated by
the same sufficiently small angle, including the corresponding
radial chain. Also a computable \(c>0\) gives the global bound
\[
                         J(z)\ge c\,e^{-|z|\sqrt{1+|z|}} .
                                                               \tag{2}
\]

These are effective estimates for the elementary class.
They apply to analytic remainder amplitudes on finite disks.
In particular they do not need an analytic continuation of the
amplitude to infinity.

## The absolute integrals and the joint lower bound away from zero

Let \(R=|z|\). First consider \(|\arg z|\le3\pi/4\) and
let \(s=\sqrt z\) with positive real part.
Then \(\Re s\ge\sqrt R/4\).
The vertical line \(v=s+it\) has exactly
\[
          p_z(s+it)=-2z^{3/2}/3-s t^2-i t^3/3 .
\]
Its absolute integral is at most a computable multiple of
\(R^{-1/4}|e^{-2z^{3/2}/3}|\).
Truncate the line at its intersections with \(|v|=H\) and
join the ends by circular arcs to the standard rays.
This is a valid deformation inside the analytic disk.
If \(x=\Re s\), the endpoint on an upper or lower
circular arc makes angle at least \(\arcsin(x/H)\) with the
imaginary axis. Throughout the joining arc the cubic real part
is at most \(-xH^2/2\). In fact in these arcs
\(\cos(3\theta)=-\sin(3u)\) for the absolute offset
\(u\in[\arcsin(x/H),\pi/6]\) from the imaginary axis,
and \(\sin(3u)\ge6u/\pi\).
The linear term has modulus at most \(RH\).
If \(\sqrt R/H\) is at most a sufficiently small absolute
constant, this is absorbed by \(xH^2/8\).
The arcs then have fixed exponential loss relative to the
saddle exponential, whose inverse modulus is bounded by
\(e^{2R^{3/2}/3}\).
Their lengths are \(O(H)\). Dividing their bound by
\(R^{-1/4}|e^{-2z^{3/2}/3}|\) costs at most a constant
times \(H R^{1/4}e^{-c\sqrt R H^2}\), which is bounded by
the elementary maximum of \(u^{1/2} e^{-cu}\) for \(u\ge0\).
All these estimates are effective.

For large \(R\) the leading Airy formula follows on this very
line. The joining arcs just estimated vanish as \(H\to\infty\)
at fixed \(z\). Expanding \(e^{-it^3/3}\) against
\(e^{-s t^2}\) gives
\[
 \begin{split}
  \operatorname{Ai}(z)&=
    (2\sqrt\pi)^{-1} z^{-1/4}e^{-2z^{3/2}/3}
                       (1+O(R^{-3/2})),\\
  \operatorname{Ai}'(z)&=
   -(2\sqrt\pi)^{-1} z^{1/4}e^{-2z^{3/2}/3}
                       (1+O(R^{-3/2})).
 \end{split}
\]
For the first estimate the linear cubic term is odd, and
\(|e^{-it^3/3}-1+it^3/3|\le t^6/18\).
For the second, integrate the factor \(-(s+it)\);
its odd \(it\) term vanishes against the unperturbed
Gaussian, and
\(|t(e^{-it^3/3}-1)|\le |t|^4/3\).
The resulting absolute Gaussian moments and
\(\Re s\ge\sqrt R/4\) give computable constants.
Hence in this sector for all effectively large \(R\),
\(J(z)\ge c R^{-1/4}|e^{-2z^{3/2}/3}|\), proving (1).

If \(3\pi/4\le|\arg z|\le\pi\), decompose the truncated chain
into minus the other two rotated chains, up to the irrelevant
choice of orientation, by adding the three paths through zero.
This is an exact chain identity at finite \(H\).
Write \(\omega=e^{2\pi i/3}\). Both rotated fields
\(z_1=\omega z,\ z_2=\omega^2 z\) now have principal arguments
of magnitude at most \(3\pi/4\).
Apply the preceding absolute-integral bounds to these two
terms. Their sum is bounded by a multiple of
\[
 R^{-1/4}\bigl(|e^{-2z_1^{3/2}/3}|+|e^{-2z_2^{3/2}/3}|\bigr).
\]
It is bounded by \(C J(z)\) for large \(R\).
Indeed in the exact connection formulas
\[
 \operatorname{Ai}(z)=-\omega\operatorname{Ai}(z_1)
                     -\omega^2\operatorname{Ai}(z_2),\qquad
 \operatorname{Ai}'(z)=-\omega^2\operatorname{Ai}'(z_1)
                     -\omega\operatorname{Ai}'(z_2)
\]
the two normalized leading derivative-to-value ratios are
opposite, each of modulus one after division by \(\sqrt R\).
This follows from
\(\omega\sqrt{\omega z}=-\omega^2\sqrt{\omega^2z}\)
on each of the two side sectors with the indicated square roots.
The ideal two-column matrix therefore has a bounded inverse
independent of \(z\).
The relative errors \(O(R^{-3/2})\) in the displayed
asymptotics are absorbed for effectively large \(R\).
The use of \((1+R)^{-1/2}\) in place of \(R^{-1/2}\)
changes the norm by bounded absolute factors.

## Bounded fields and finite endpoint changes

For (2), apply the Airy differential equation along the segment
from zero to \(z\) to the pair
\((\operatorname{Ai},\operatorname{Ai}'/\sqrt{1+|z|})\).
The matrix has induced sum norm at most \(\sqrt{1+|z|}\).
Gronwall's inequality in the reverse direction gives (2)
with \(c=\operatorname{Ai}(0)>0\).
This also bounds \(J(z)\) uniformly below on the remaining
compact \(R\)-disk. The two standard rays bound the integral
there directly by the integral of
\(2\exp(-v^3/3+R v)\) for real \(v\ge0\).
Thus all constants needed for (1) on the bounded fields
are computable.

For a simultaneous small rotation of the truncated endpoints,
join each endpoint to the unrotated ray on its circle.
Both arcs stay strictly inside their cubic decay sectors.
They contribute at most \(C H e^{-c_2 H^3}\) if
\(|z|\le\delta H^2\) and \(\delta\) is decreased as needed.
Bound (2) absorbs this term into (1) on a further fixed
decrease of \(\delta\).
The same comparison absorbs any fixed finite exterior
contour with phase loss \(e^{-c_3 j}\) in a local
Airy calculation with \(H\) proportional to \(j^{1/3}\),
provided its scaled field disk is small enough.
The absolute constants in these decreases follow from (2).

The Airy contour normalization and sign of \(\operatorname{Ai}(0)\)
are the standard real normalization.
[Airy contour formula](https://dlmf.nist.gov/9.5.E4)
