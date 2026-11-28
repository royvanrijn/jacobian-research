# An explicit initial interval for every higher cell count

Use the equation labels of the [interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md).
Put \(p_k=\gamma_H(F_{12})\), where \(H\) is the reference pair
plane and \(\gamma_H\) is standard Gaussian probability measure
on \(H\). Thus \(v_k=\phi_0p_k\) and \(\ell=k\phi_0p_k\). For every
\(k\ge6\), \(d\ge k-1\) and
\[
      0<\rho\le\frac{p_k}{2-p_k},
\]
the full balanced interface form has a positive dimension-independent
gap modulo rotations. This interval contains
\(0<\rho\le 1/[k(k-1)]\) by the facet-mass estimate in (S81).

Use the normalization in (FV2) of the
[five-cell standard argument](FIVE_CELL_STANDARD_GAP.md), replacing
\(5\) by \(k\) and using the \(k-2\) adjacent terms. In the active
space, write \(\sigma=\sqrt{1-\rho^2}\), \(\mu=\gamma_H|_F/p_k\),
\(D_\rho=\sigma a_\rho/(\rho\ell)\ge1\) by (B7), and let \(A_\rho\)
denote the standard interaction in this normalization. For an ordinary
balanced invariant \(h\), set \(w=(p_0/|p_0|)\cdot x\) and
\(m=\int hw\,d\mu\). By (S77), the rank-one degree-one
interaction term of \(A_\rho\) is \(\frac12|m|^2\).
Also \(|m|^2\le\|h\|_\mu^2\) by the same convex covariance bound
used in (S77), since \(m=\int h(w-\int w\,d\mu)\,d\mu\).

For the ambient Hermite degrees at least two, an adjacent term
is bounded above by the same-facet term degree by degree, by
Cauchy–Schwarz and orthogonality of the transposition within
each chaos. The sum of the kernel weights equals \(k\), so the
standard normalization leaves a same-facet coefficient
\(\sigma/\rho\). Extend \(h\) by zero onto \(H\). A planar Hermite
coefficient of degree \(j\ge2\) has multiplier \(\rho^{j-1}\) after
summing the even normal terms, using
\(\sum_{r\ge0}\binom{2r}{r}\rho^{2r}/4^r=\sigma^{-1}\).
For \(j=1\) the ambient degree-one contribution has already been
removed, leaving multiplier \(1-\sigma\). The planar degree-zero
coefficient vanishes by ordinary balance. Thus \(\rho\) is a
cap on all the tail multipliers, with a correction
\(\rho-(1-\sigma)\) on the linear coefficient. The zero extension
has squared norm \(p_k\|h\|_\mu^2\) and its planar linear coefficient
has square \(p_k^2|m|^2\), because the moment vector lies on the
axis \(p_0\). The interaction is integrated against \(\mu^2\),
so its two terms carry the factor \(p_k^{-2}\). Adding the
degree-one interaction back, we obtain
\[
 \langle h,\sigma A_\rho h\rangle_\mu
 \le \frac{\rho}{p_k}\|h\|_\mu^2
             +(1-\rho-\sigma/2)|m|^2.
\]
For \(0<\rho<3/5\), the second coefficient is positive. Dividing
the standard form by the facet norm yields a gap at least
\[
    \frac{\ell}{\sigma}
       \left[\frac{\sigma}{2}-\rho\left(\frac{1}{p_k}-1\right)\right].
\]
The bracket decreases in \(\rho\). At the specified endpoint it equals
\[
       \frac{\sqrt{1-p_k}-(1-p_k)}{2-p_k}>0 .
\]

For the other blocks, (S56) and (S87) already handle rotations
and all nonstandard two-row shapes. Every other representation
has Hermite floor \(D_\lambda\ge2\) by (S80), because a shape
with at least three rows and floor below two must be the
rotation hook. Use \(b=\sqrt{2/3}\) in (S81) to bound each of
these gaps below by \(\ell/\sigma-6\sqrt3\,\phi_0\rho\).
The facet probability is the probability that a \(N(0,1/2)\)
variable exceeds \(k-2\) independent standard normals, upon
removing the common trace. It therefore decreases with \(k\).
Formula (FV1) gives \(p_k\le p_5<1/4\) since
\(\arcsin(1/3)<\pi/6\). The endpoint is below \(1/7\), hence
within \(b\) and \(3/5\), and
\[
    6\sqrt3\,\phi_0\rho
       \le \ell\frac{6\sqrt3}{k(2-p_k)}
       <\ell\frac{4\sqrt3}{7}<\ell .
\]
Taking the minimum of the gaps over the active blocks and applying
(S42) gives the assertion in all ambient dimensions.
