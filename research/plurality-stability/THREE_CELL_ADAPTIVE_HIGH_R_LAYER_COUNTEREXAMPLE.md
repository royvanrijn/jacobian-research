# Sharp adaptive energy near one need not control tagged perimeter

The sharp high-correlation limit for the adaptive
surface energy does not imply strict perimeter
convergence of its near-minimizers. There are ordinary
balanced partitions converging to the centered simplex
in Gaussian label distance, with asymptotically minimal
unrestricted adaptive energy near correlation one,
whose Gaussian perimeters have a strictly larger
limit than the simplex perimeter. The partitions
have a thin, rapidly alternating layer of the third
label on a fixed segment of a simplex wall.

[Adaptive high-correlation lower bound and minimum](THREE_CELL_ADAPTIVE_SURFACE_HIGH_CORRELATION_LIMIT.md)
[The unrestricted adaptive functional](THREE_CELL_ADAPTIVE_SURFACE_UNRESTRICTED_DIRECTIONS.md)

## Pointwise direction optimization

For an ordinary partition let
\(H_e=T_s\nu_e\), where \(s=\sqrt r\) and each
\(\nu_e\) is the Gaussian tagged interface measure.
At every observation \(x\), the contribution of
directions on that tag is a Hilbert vector of norm
at most \(H_e(x)\). Thus
\[
 F_\infty(A)\le
       \|{\cal V}(H_1,H_2,H_3)\|_{L^2(\gamma_2)},\qquad
 {\cal V}(h)=\max_{|v_e|\le h_e}
                       \Big|\sum_e d_e\otimes v_e\Big| .     \tag{HA1}
\]
Use cyclically oriented label vectors
\(d_1=e_1-e_2,d_2=e_2-e_3,d_3=e_3-e_1\).
They have length \(\sqrt2\) and pairwise inner product
\(-1\), so the squared norm inside the maximum is
\(3\sum_e |v_e|^2-|\sum_e v_e|^2\).
Keeping the other vectors fixed, convexity makes the
maximum on each directional ball attainable on its
boundary. The smallest norm of the sum of vectors
with those three lengths is the excess of the largest
one over the sum of the other two, or zero when there
is no excess. In particular for \(a\ge2b\ge0\),
\[
                    {\cal V}(a,b,b)^2=2(a+b)^2.              \tag{HA2}
\]
The formula also holds when the dominant tag is a
different pair. The map \({\cal V}\) is Lipschitz
with constant \(\sqrt2\) for the sum of the absolute
capacity differences, by replacing each direction
by a vector with the new maximal length on the same
ray.

## Inserting a thin layer on a fixed wall segment

Let \(T\) be the centered regular simplex. Fix an
interior straight segment \(I\) of one ray,
separated by a positive distance from the other
two rays and the vertex. Work in orthonormal
coordinates \(z,n\) on the tangent and normal axes
of the chosen ray, so \(I=[a,b]\times\{0\}\) with
\(0<a<b\). Write \(i,j\) for its two labels and \(k\)
for the spectator. Fix \(0<\theta<1/3\).

Put \(\ell=\sqrt{1-r}\) and
\(\sigma=\sqrt{1-r^2}\). Divide \([a,b]\) into
equal intervals of length \(h=h_r\) of order
\(\ell^2\), and let \(J_h\) be the union of their
first \(\theta\) fractions. In \(J_h\times(-\delta,\delta)\),
with \(\delta=\ell^4\), assign the spectator \(k\)
in place of the adjacent labels. Thus the original
ray segment on \(J_h\) is replaced by the pair
walls at \(n=\pm\delta\), and by short cap segments
at the ends of the small intervals. The caps
have total length \(O(\delta/h)\). The two
independent volume errors are \(O(\delta)\).
Correct them exactly by smooth compact flows on
fixed patches of other regular simplex walls,
disjoint from a neighborhood of \(I\). The
Gaussian flux derivative matrix is invertible
because the chosen patches transfer between two
independent label pairs. The implicit function
theorem gives flow times \(O(\delta)\).
Call the resulting ordinary balanced partition
\(A_r\). Its label distance to \(T\) tends to
zero. Its perimeter satisfies
\[
      P_\gamma(A_r)\ \longrightarrow\
            P_\gamma(T)+\theta\int_I\phi_2\,d\mathcal H^1 .
                                                               \tag{HA3}
\]
Indeed the removed Gaussian length tends to
\(\theta\int_I\phi_2\), each new parallel segment
has that same limit, and the cap and balancing
costs tend to zero.

## Energy of the alternating layer

The smoothed tagged arclength at observation \(x\)
is the integral of the Lebesgue conditional
Gaussian density with mean \(sx\) and covariance
\(\ell^2 I\) over that interface. Put
\(g_\ell(t)=(2\pi\ell^2)^{-1/2}\exp(-t^2/(2\ell^2))\).
For the whole chosen reference ray its smoothed
arclength is
\[
 P(x)=g_\ell(sn)\,
                   (g_\ell*1_{[0,\infty)})(sz).
\]
For the segment \(I\) the same formula with
\(1_{[a,b]}\) gives \(Q(x)\), where \(0\le Q\le P\).
The layer construction has, to leading order
in \(L^2(\gamma_2)\) multiplied by \(\sqrt\sigma\),
the following contribution from the chosen ray:
\[
                 (P-\theta Q,\ \theta Q,\ \theta Q)
                                                               \tag{HA4}
\]
on tags \(ij,ik,kj\), respectively. In fact the
antiderivative of \(1_{J_h}-\theta1_{[a,b]}\)
is supported in \([a,b]\) and bounded in
absolute value by \(h\). Integration by parts
therefore bounds the uniform tangential
convolution error by \(C h/\ell\). The factor
\(\|g_\ell(sn)\|_{L^2(\gamma_1)}\) is \(O(\ell^{-1/2})\).
The displaced normal kernels have
\(L^2(\gamma_1)\) difference from \(g_\ell(sn)\)
bounded by \(C\delta\ell^{-3/2}\), by the
derivative bound and translation integral. These
estimates give the claimed scaled error.

The caps have scaled error tending to zero as
well. For a compact source point \(y\) the
\(L^2(\gamma_2)\) norm of the point-source
Mehler kernel is \(M_r(y,y)^{1/2}=O(\ell^{-1})\).
Thus the cap contribution multiplied by
\(\sqrt\sigma\) is \(O(\delta/(h\sqrt\ell))=o(1)\).
For the compact exterior balancing flows, the
same norm bound and its differentiated kernel
bound \(O(\ell^{-2})\) give a scaled norm change
\(O(\delta\ell^{-3/2})=o(1)\). The derivative
bound follows by differentiating
\(\langle M_s(\cdot,y),M_s(\cdot,y')\rangle=M_r(y,y')\)
in both kernel arguments on the diagonal.
The Gaussian surface Jacobians of the compact
flows change by \(O(\delta)\), which is also
covered by this estimate. The unmodified rays
retain their original smoothed arclengths.
Thus (HA4) and these estimates approximate the
entire capacity vector in the scaled norm.

In a fixed small tubular neighborhood of \(I\),
the other two reference rays have negligible
scaled smoothed arclength, since they stay a
fixed positive distance from that neighborhood.
Outside it the segment contribution \(Q\)
has negligible scaled norm by the same Gaussian
distance bound. In the neighborhood,
\(P-\theta Q\ge 2\theta Q\), so by (HA2) the
pointwise optimized norm for (HA4) is exactly
\(\sqrt2 P\). It is the same as that of the
unmodified ray alone. The Lipschitz bound of
(HA1) then gives
\[
 \sqrt{\sigma}F_\infty(A_r)
   \ \le\
    \sqrt{\sigma}\|{\cal V}(H(T))\|_2+o(1),
                                                               \tag{HA5}
\]
where \(H(T)\) is the reference triple of smoothed
ray measures. For the unmodified rays their
mixed products have negligible integrals after
multiplication by \(\sigma\), since distinct rays
only meet at the vertex and their joint Mehler
cross integrals are \(O(1)\). The diagonal terms
satisfy the usual one-ray limit. Therefore
\[
 \limsup_{r\uparrow1}\sigma F_\infty(A_r)^2
                          \le \sqrt{2/\pi}\,P_\gamma(T).
\]
The sharp high-correlation lower bound gives
the reverse inequality, since \(A_r\to T\) in
Gaussian label measure. This proves asymptotic
minimality despite the larger perimeter in
(HA3). In particular label compactness together
with the sharp energy limit leaves room for
collapsed tagged interfaces in highly oscillatory
near-minimizing sequences.
[The simplex ray integrals and the sharp minimum](THREE_CELL_ADAPTIVE_SURFACE_HIGH_CORRELATION_LIMIT.md)
