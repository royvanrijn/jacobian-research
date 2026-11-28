# Conformal variation of Gaussian conical partitions

On the unit sphere \(S^{d-1}\), use the angular noise operator
\(P_r\) induced by the directions of a standard Gaussian pair at
scalar correlation \(0<r<1\). For \(d\ge3\), every conical noise
stability maximizer with at least two prescribed positive cell
masses has a nonzero cell centroid. The proof uses a strict positive
trace of its priced conformal Hessian before the volume constraint.
Centroid cancellation does not, by itself, give a sign for the
individual directions of this form, even on ordinary balanced
three-cell partitions. This bounds what the conformal calculation
alone gives toward same-dimension reduction of conical optimizers.

## Angular operator and conformal derivative

Write \(d\omega\) for uniform probability measure on the sphere,
and let \(u=(u_i)\) be the angular indicators of the conical cells.
Their stability is \(\sum_i\langle u_i,P_r u_i\rangle\). The kernel
of \(P_r\) is a smooth function of the dot product of its two
directions. Indeed it is obtained from the joint Gaussian kernel
by integrating the two positive radii with their polar
Jacobians. The covariance quadratic form has a uniform positive
lower bound in the radial variables, which permits differentiation
of the angular kernel to every order. In particular all the
distributional pairings below are well defined. This is the
angular operator in the spherical Gaussian-noise decomposition.
[Hwang, Neeman, Parekh, Thompson and Wright, Sections 4 and 5](https://people.eecs.berkeley.edu/~jswright/papers/quantum-borell.pdf).

Write \(\lambda_\ell\) for its eigenvalue on harmonics of degree
\(\ell\), with \(\lambda_0=1\), and \(L_\ell=\ell(\ell+d-2)\)
for the corresponding eigenvalue of the negative spherical
Laplacian. For the Gaussian angular operator there is a common
random clock \(\tau\), positive and finite almost surely, such that
\[
                \lambda_\ell=\mathbb E e^{-L_\ell\tau}.       \tag{CF1}
\]
To see this, let \(Y=r(X+B_{r^{-2}-1})\) where \(X\) is a
standard Gaussian and \(B\) is an independent standard Euclidean
Brownian motion. Conditional on the starting radius of \(X\), the
direction of \(X+B_t\) is spherical Brownian motion, with generator
\(\Delta_S/2\), at an independent time change determined by the Bessel
radius. Averaging over the initial Gaussian radius preserves
independence from the angular Brownian motion conditional on the
initial direction. Take \(\tau=\tfrac12\int_0^{r^{-2}-1}R_s^{-2}ds\)
for this Bessel radius \(R_s\), which starts at a positive radius
and does not hit zero in dimension \(d\ge3\). The final positive
scale factor \(r\) does not affect direction. This proves (CF1).
[Byczkowski, Chorowski, Graczyk and Małecki, Section 6.3, for the radial–angular Brownian skew product](https://arxiv.org/pdf/1109.0182).

For a constant vector \(v\), put \(l_v(\omega)=v\cdot\omega\)
and let \(D_v=X_v\cdot\nabla_S\), with
\(X_v(\omega)=v-l_v(\omega)\omega\). This vector field generates
a smooth flow on the sphere. Its divergence is
\(-(d-1)l_v\). Thus the volume derivative of the \(i\)-th angular
cell under its flow is \(-(d-1)\int l_v u_i\,d\omega\). The
Euclidean Gaussian cell centroid is a positive multiple of
\(\int\omega u_i\,d\omega\), by radial independence.

Suppose the tuple is an angular winner for the scores
\(s_i=P_r u_i-\beta_i\) with constant prices. These are the
stationary scores at a conical maximizer by the angular
linear assignment argument of (CH6). Differentiating the
pointwise winner identity gives \(\sum_i s_i D_v u_i=0\)
distributionally. Let \(u_i^t\) be the transported indicator
under the conformal flow, and subtract twice the priced mass
from its stability. The resulting second variation divided
by two at zero is
\[
 \mathcal C(v)=\sum_i\bigl(\langle D_v u_i,P_r D_v u_i\rangle
                    -\langle D_v u_i,D_vP_r u_i\rangle\bigr).
                                                               \tag{CF2}
\]
Indeed \(u_i^t\) has first distributional derivative \(-D_vu_i\)
and second derivative \(D_v^2 u_i\). Gaussian radial averaging
has already been absorbed into the angular kernel. In sphere
integration by parts, \(D_v^*=-D_v+(d-1)l_v\), and the last
multiplication term cancels in the sum by the winner identity.
Differentiation can equally be done by changing variables in
the smooth angular kernel. If the centroid projections onto
\(v\) vanish, the volumes can be repaired at arbitrarily small
quadratic cost by the angular tree transfers in (CH7). The
proof uses only the stationary angular scores and the continuity
of the transported scores. Consequently actual conical maximality
requires \(\mathcal C(v)\le0\) for every such \(v\).
[Angular score rule and exact angular repairs](THREE_CELL_CONICAL_HOMOGENIZATION.md#equality-of-consecutive-conical-optima).

## Trace and exclusion of zero centroids

For \(\ell\ge1\), define \(\delta_\ell=\lambda_\ell-\lambda_{\ell+1}\)
and also \(\delta_0=1-\lambda_1\). From (CF1), the decreasing
function \(x\mapsto \tau e^{-x\tau}\) gives
\[
 \frac{\delta_{\ell-1}}{2\ell+d-3}
              >\frac{\delta_\ell}{2\ell+d-1}                 \tag{CF3}
\]
by averaging over adjacent intervals of the Laplacian spectrum.

Let \(u_{i,\ell}\) be the degree-\(\ell\) spherical harmonic part.
For any orthonormal ambient basis \((e_j)_{j=1}^d\) the trace is
\[
 \sum_j\mathcal C(e_j)=
 \sum_{i,\ell\ge1}\frac{\ell(\ell+d-2)}{2\ell+d-2}
   \bigl[(\ell+d-2)\delta_{\ell-1}-\ell\delta_\ell\bigr]
                                     \|u_{i,\ell}\|_2^2 .      \tag{CF4}
\]
The calculation is as follows. For a degree-\(\ell\) homogeneous
harmonic polynomial \(p\) restricted to the sphere, multiplication
by \(l_v\) decomposes as \(L_vp+U_vp\), where
\(L_vp=(\partial_vp)/(2\ell+d-2)\) has degree \(\ell-1\) and
\(U_vp\) has degree \(\ell+1\). In Euclidean coordinates the
upper polynomial is \((v\cdot x)p-|x|^2L_vp\), which is harmonic
by direct differentiation. Tangential differentiation is
\(D_vp=(\ell+d-2)L_vp-\ell U_vp\). Moreover
\[
 \sum_j\|L_{e_j}p\|_2^2=\frac{\ell}{2\ell+d-2}\|p\|_2^2,\qquad
 \sum_j\|U_{e_j}p\|_2^2=\frac{\ell+d-2}{2\ell+d-2}\|p\|_2^2.
\]
The first equality uses \(\sum_j\|\partial_{e_j}p\|_2^2
=(\ell^2+L_\ell)\|p\|_2^2\), by the radial and tangential
decomposition of the gradient. The second uses orthogonality
and \(\sum_j l_{e_j}^2=1\). The possible mixed terms after
summing over \(j\) pair \(L_{e_j}p_{\ell+1}\) with
\(U_{e_j}q_{\ell-1}\). The former is orthogonal to the lower
part of \(l_{e_j}q_{\ell-1}\). Pairing it with the full product
and summing over \(j\) leaves a multiple of
\(\langle p_{\ell+1},q_{\ell-1}\rangle=0\) by Euler's identity.
Equation (CF2) is \(\sum_i\langle D_v u_i,[P_r,D_v]u_i\rangle\).
Insert these lower and upper terms to obtain (CF4). One may
first calculate on finite harmonic expansions and pass to the
limit, since after integration by parts (CF2) uses only smooth
kernels and is continuous on \(L^2\).

For \(d\ge3\) the coefficients in (CF4) are strictly positive
for every \(\ell\ge1\). Indeed (CF3) suffices since
\[
 (\ell+d-2)(2\ell+d-3)-\ell(2\ell+d-1)
                         =(d-3)(2\ell+d-2)\ge0 .
\]
A partition with positive masses has nonconstant indicators,
so the trace is positive. If all centroids were zero, every
direction would be volume tangent. Angular maximality would
give \(\mathcal C(e_j)\le0\) for all \(j\), a contradiction.

## Obstruction to an individual directional sign

The obstruction concerns the quadratic expression (CF2),
extended algebraically to arbitrary indicator or fractional
tuples even when there is no stationary winner rule. In any
dimension \(d\ge3\), take a nonzero real harmonic polynomial
of even degree \(\ell\ge2\) independent of \(x_1\), for example
\(p(x)=\operatorname{Re}(x_2+i x_3)^\ell\), with \(v=e_1\).
Then \(L_v p=0\). The degree of \(D_vp=-\ell l_v p\) is
purely \(\ell+1\), so
\[
    \langle D_vp,[P_r,D_v]p\rangle
          =-\delta_\ell\|D_vp\|_2^2<0.                      \tag{CF5}
\]
For small positive \(\epsilon\), this gives a fractional
balanced tuple \((1/3+\epsilon p,1/3-\epsilon p,1/3)\)
with negative (CF2) along \(v\) and no centroids. It can also
be approximated weakly by ordinary balanced zero-centroid
tuples. Partition the antipodal quotient of the sphere into
finer patches with diameters tending to zero, and split each
patch by the three exact masses obtained by integrating the
fractional tuple on it. The quotient measure is nonatomic.
The lifted indicators are even, which gives zero centroids,
and converge weakly to the fractional tuple by testing on
continuous functions. The distributional expression (CF2)
is weakly continuous on bounded tuples because after moving
the derivatives to the kernel it is a compact quadratic
form with a smooth kernel. Thus some ordinary balanced even
tuples also give strict negativity of that expression.

These latter tuples need not be stationary, so (CF5) does not
claim the existence of a stable nonsimplex cone. It shows that
the centroid condition and indicator constraints do not alone
give a positive directional sign by spectral diagonalization.
The remaining rank reduction needs a stronger use of the
angular winner equation or a different volume-preserving
test. In particular the full trace in (CF4) includes the
directions that change masses whenever the centroid span is
nonzero.
