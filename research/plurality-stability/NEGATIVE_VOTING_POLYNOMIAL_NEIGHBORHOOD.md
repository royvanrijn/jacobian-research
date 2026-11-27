# A polynomial local neighborhood for negative Gaussian voting

Let \(P_r(F,G)=\sum_i\langle F_i,T_{-r}G_i\rangle\) for two
measurable fractional Gaussian assignments with corresponding
masses equal. The common mass vector is free. Let \(d_T(F,G)\)
be the infimum of their summed label \(L^1\) distance to a common
cylindrical regular simplex.

There are absolute \(A>0\) and \(k_0\) such that for every \(k\ge k_0\),
\(d\ge k-1\), \(0<r\le1/(k-1)\), and \(d_T(F,G)\le k^{-A}\),
\[
 P_r(F,G)\ge P_r(T,T).
 \tag{PN1}
\]
Equality requires both assignments to be the same simplex, up to
simultaneous orthogonal maps, labels and null sets. The proof uses
[the weighted voting-scale Hessian gap](NEGATIVE_VOTING_HESSIAN.md#the-full-negative-paired-gap)
and projects measurable transfers to the facets with polynomial
error bounds. It preserves the exact mass interaction.

## Global score and centered kernel bounds

Write \(E=\mathbf1^\perp\), \(n=k-1\), and split the ambient
coordinates as \((x,z)\in E\oplus\mathbb R^{d-n}\).
Set \(n_{ij}=(e_i-e_j)/\sqrt2\), \(H_{ij}=n_{ij}^\perp\cap E\),
and
\[
 V_{ij}=-T_{-r}(1_{T_i}-1_{T_j}),\qquad
 a_r(y)=\partial_{n_{ij}}V_{ij}(y),\quad
 A_r(y)=a_r(y)/(r\ell)\qquad (y\in H_{ij}),
\]
where \(\ell=k v_k\) for active facet area \(v_k\).
Throughout, a polynomial factor is at most \(C k^B\) with
absolute \(C,B\) independent of \(d,r\).

For \(x=y+t n_{ij}\), \(t\ge0\), the reflection integral
[(N4)](NEGATIVE_SMALL_CORRELATION.md#exact-negative-correlation-scores)
and \(\sinh u\ge u\) for \(u\ge0\) give
\[
 V_{ij}(y+t n_{ij})\ge t e^{-\beta t^2}a_r(y),
 \qquad \beta=\frac{r^2}{2(1-r^2)}.
 \tag{PN2}
\]
Indeed, the conditional normal coordinate on \(T_i\) is nonnegative.
In the reflected difference integral the normal factor is a positive
hyperbolic sine times \(e^{-\beta t^2}\). The linear term of
that sine gives the derivative at zero.

The row formula used in
[(VH5)](NEGATIVE_VOTING_HESSIAN.md#the-reciprocal-local-weight)
applies on the entire pair plane. The distributional normal
derivative of \(1_{T_i}-1_{T_j}\) has coefficient 2 on the \(ij\)
facet and \(1/2\) on each \(il\) and \(jl\) facet for \(l\notin\{i,j\}\).
At each \(y\in H_{ij}\), reflection identifies the two latter
kernel integrals. The same Jensen estimate as on the facet yields
\[
 A_r(y)\ge \sigma^{-1}
       \exp\!\left(-\beta|y|^2-\frac{r m_k}{2\sigma}p_{ij}\cdot y\right)
       \ge c e^{-C|y|^2/k},
 \qquad \sigma=\sqrt{1-r^2},
 \tag{PN3}
\]
where \(p_{ij}\) is the unit tie axis and \(m_k=O(\sqrt{\log k})\)
by (VH6). Absorb the linear term by Young's inequality, using
\(k r^2 m_k^2=O(1)\). Differentiation in one unit direction of the
Gaussian convolution for a bounded indicator difference also gives
\(|\partial_{n_{ij}}V_{ij}|\le Cr\). Thus \(A_r\le C/\ell\) and
\(|V_{ij}(y+t n_{ij})|\le Cr |t|\). In particular,
\[
 U_{ij}(x):=V_{ij}(x)/(r\ell)
       \ge c e^{-C|x|^2/k}\min(t,1),\qquad t\ge0 .
 \tag{PN4}
\]

For the active Mehler density relative to product Gaussian measure,
differentiation in correlation gives
\[
 |M^x_{-r}(x,y)-1|\le r\,\mathrm{poly}(k)\,W(x)W(y),
 \qquad W(x)=(1+|x|^2)e^{C_1|x|^2/k},
 \tag{PN5}
\]
with absolute \(C_1\) sufficiently large. The Mehler prefactor is
bounded on the voting interval, its exponent is bounded above by
\(C(|x|^2+|y|^2)/k\), and its logarithmic derivative in correlation is
bounded by a polynomial in \(k,|x|,|y|\) of degree at most two in the
spatial variables. On \(|x|,|y|\le K\sqrt{k}\), the spatial derivatives
of the same kernel are bounded by \(r\,\mathrm{poly}(k)\) for fixed \(K\).

Let \(P_0\) and its factor versions be Gaussian expectations.
The exact decomposition of the centered operator is
\[
 T_{-r}-P_0=(T^x_{-r}-P^x_0)\otimes T^z_{-r}
                 +P^x_0\otimes(T^z_{-r}-P^z_0).
 \tag{PN6}
\]
The bilinear norm of \(T^z_{-r}-P^z_0\) on
\(L^{3/2}\times L^{3/2}\) is at most \(2r\).
Factor \(T^z_{1/\sqrt2}\) on both sides of \(T^z_{-2r}-P^z_0\),
then use the Gaussian \(L^{3/2}\)-to-\(L^2\) contraction and the
middle \(L^2\) norm \(2r\). The same argument without centering
gives norm at most 1 for \(T^z_{-r}\).
[Nelson's hypercontractivity as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573).
By (PN5)–(PN6) the absolute centered interaction of two signed
transfers is at most \(r\,\mathrm{poly}(k)\) times the product of
the \(L^{3/2}(\gamma_{d-n})\) norms of their \(W\)-weighted
absolute active slice masses.

## Active-slice bounds and normal projection

For an ordered transfer \(h=1_{T_i}F_j\) or \(1_{T_i}G_j\),
\(i\ne j\), let \(c_h(z)=\int_E U_{ij}h\,d\gamma_n\) and
\(j_h(z)=\int_E W h\,d\gamma_n\). Weighted Cauchy–Schwarz gives
\[
 j_h(z)^2\le \mathrm{poly}(k)c_h(z).
 \tag{PN7}
\]
For \(t\ge1\) this follows directly from (PN4).
For \(0\le t<1\), disintegrate along \(n_{ij}\) and use
\((\int h\phi(t)dt)^2\le2\phi(0)\int t h\phi(t)dt\). This is
the minimum first-moment bound for a density bounded by
\(\phi(0)\) on the positive axis. Next use tangential
Cauchy–Schwarz with the reciprocal of the score weight.
The required integrals of \(W^2(x)e^{C|x|^2/k}\), in the active
space or on the pair plane, are polynomially bounded by Gaussian
radial integration in dimension at most \(k\).
With the additional restriction \(|x|>R=K\sqrt{k}\), for a
sufficiently large absolute \(K\), the right side of (PN7)
can be multiplied by \(e^{-c k}\). For \(t<1\) the plane
radius then exceeds \(\sqrt{R^2-1}\), and the same
Gaussian radial tail estimate applies.

Choose the simplex whose moment matrix is closest to the average
moment matrix of \(F,G\), and use its active coordinates. The
average moment error has zero projection onto the tangent of the
regular moment orbit. If the original orbit distance is \(\eta>0\),
the distance \(\eta_*\) to the aligned simplex is at most
\(\mathrm{poly}(k)\sqrt{\eta}\). Indeed an initially aligning
simplex at distance at most \(2\eta\) has average moment error
\(O(\sqrt\eta)\) by first Hermite projection. The nearest moment
simplex differs in its moment matrix by \(O(\sqrt\eta)\) from
the first reference. Their corresponding unit pair normals differ
by at most \(C\sqrt\eta/\ell\), and the union bound on Gaussian
wedges between pair halfspaces gives the claim. Distance zero
already gives a common simplex by compactness of its orbit in
each fixed ambient dimension.

Set \(J=\sum_h\int c_h\,d\gamma_{d-n}\) across both assignments.
The linear part of the exact deficit is \(r\ell J\).
If \(J=0\), (PN4) forces both assignments to be the reference
simplex. Otherwise normalize by \(\sqrt J\). The bound
\(|V_{ij}|\le Crt\), with \(t\le |x|\) and
\(\sum_{i\ne j}h\le1\) on each assignment, gives
\(J\le \mathrm{poly}(k)\sqrt{\eta_*}\).
Thus \(J\) can be made polynomially small of any fixed order
by taking \(\eta\le k^{-A}\) with \(A\) sufficiently large.

Retain only transfers with \(|x|\le R\), \(0\le t\le\xi\le1\)
and \(y=x-t n_{ij}\in F_{ij}\), and project normally to \(y\),
preserving \(z\). The discarded \(W\)-weighted active slice
masses, in \(L^{3/2}\) after division by \(\sqrt J\), are bounded by
\[
 \mathrm{poly}(k)
   \left(e^{-c k}+\sqrt{\xi}+(\sqrt J/\xi)^{1/3}\right).
 \tag{PN8}
\]
The exterior estimate in (PN7) accounts for the first term.
For \(|x|\le R,\ t>\xi\), (PN4) gives an \(L^1\) slice mass
bound \(\mathrm{poly}(k)J/\xi\). Interpolation with the \(L^2\)
bound from (PN7) gives the third term. For
\(|x|\le R,\ t\le\xi,\ y\notin F_{ij}\), some untied label \(l\)
satisfies \(0<y_l-y_i\le\xi/\sqrt2\), because \(x\in T_i\).
Each \(y_l-y_i\) on the tie plane has Gaussian variance bounded
below by an absolute constant. The union of the projected strips
therefore has plane measure at most \(Ck\xi\). Use the same
normal mass-cost inequality and tangential Cauchy–Schwarz as in
(PN7), restricted to those strips. On the cap all the weights
are polynomial, giving the middle term of (PN8) in \(L^2\)
and hence in \(L^{3/2}\). The number of ordered labels costs only
a polynomial factor.

Form two signed facet speeds \(f,g\) from the projected transfers,
by subtracting opposite one-sided densities. The one-sided
densities relative to \(d\nu_{ij}(y)d\gamma_{d-n}(z)\) are
\[
 \frac1{\phi(0)\sqrt J}\int_0^\xi h_{ij}^{\rm kept}
                  (y+t n_{ij},z)\phi(t)\,dt .
\]
With \(D_A\) the weighted squared norm with multiplier \(A_r\),
the sharp one-dimensional mass-cost inequality and (PN2) give
\[
 D_A(f,f)+D_A(g,g)\le 2e^{\beta\xi^2}.
 \tag{PN9}
\]
For each pair, the squared difference of opposite one-sided
densities is at most the sum of their squares. Their weighted
cost is bounded by \(e^{\beta\xi^2}\) times the retained exact cost.

Let \(K_r^0=K_r-\langle C\cdot,C\cdot\rangle\) be the centered
interface form, where \(C\) records the label masses of facet
speeds. Projection preserves the kept mass of each active slice,
so the second term in (PN6) is unchanged. For its first term,
normal projection changes the active kernel on the cap by
\(r\,\mathrm{poly}(k)\xi\). The slice bounds in (PN7) control
the total retained product. They also control the projected
slice masses because orthogonal projection decreases \(|x|\).
The discarded interactions are bounded by (PN6)–(PN8).
Consequently
\[
 \frac{\sum_i\langle F_i-T_i,(T_{-r}-P_0)(G_i-T_i)\rangle}
              {r\ell J}
 = \frac{K_r^0(f,g)}{r\ell}+O(\varepsilon_k),\qquad
 \varepsilon_k=\mathrm{poly}(k)
    \left(e^{-c k}+\sqrt{\xi}+(\sqrt J/\xi)^{1/3}\right).
 \tag{PN10}
\]
Dividing by \(\ell\asymp\sqrt{\log k}/k\) changes only the
polynomial factors.

## Centered Hessian and the volume mass term

Write \(m_i=\int(F_i-T_i)=\int(G_i-T_i)\).
The discarded-transfer bounds give
\(Cf,Cg=m/\sqrt J+O(\varepsilon_k)\).
Also \(\Pi_{\rm tan}\mathcal D(f+g)=O(\varepsilon_k)\), where
\(\mathcal D\) is the facet first-moment map and
\(\Pi_{\rm tan}\) projects onto the regular moment tangent space.
For active moments the weight \(W\) bounds the discarded part
and projection moves kept points by at most \(\xi\). For unused
moments projection is exact, and the discarded part is bounded
by Hölder against a unit Gaussian linear functional in \(L^3\),
then taking the supremum over such tests for each label.
The normalization uses the exact average volume moment gauge.

Restore the summed moment gauge by subtracting a common
rotational speed. The operator \(\mathcal D\) is polynomially
bounded in the \(D_A\) norm, uniformly in \(d\). For each row
test against a unit vector and use Cauchy–Schwarz with (PN3)
and the Gaussian surface integrals of the squared linear
functional. The transverse part uses only its one-dimensional
second moment. The inverse from the tangent space to rotation
speeds is also polynomially bounded. In fact the tangent
matrix of a rotation \(\Omega\) is the active regular moment
matrix composed with \(\Omega^T\). Its regular singular value
\(\ell/\sqrt2\) controls the active and mixed blocks of
\(\Omega\), which determine the speed. Integrating their
linear normal components against \(A_r\,d\nu\,d\gamma_{d-n}\)
uses only \(A_r\le C/\ell\) and Gaussian surface second
moments. This correction preserves \(Cf,Cg\).

The centered interface operator \(K_r^0/(r\ell)\) has polynomial
norm in the \(D_A\) space, uniformly in \(r,d\).
In the active space bound the first degree by weighted
Cauchy–Schwarz and (PN3), and apply
[(VH8)](NEGATIVE_VOTING_HESSIAN.md#the-interaction-beyond-degree-one)
to degree two and above. In unused degree \(j\ge1\), the
active interaction has factor \((-r)^j\), which absorbs the
factor \(r^{-1}\) of the normalized active constant term.
The latter is the squared active mass map on that unused
degree. The actual full constant term occurs only in unused
degree zero, where it was removed. This proves the claimed
bound by transverse Hermite expansion. Consequently the
rotation correction changes the centered interaction in (PN10)
and the squared norm in (PN9) by \(O(\varepsilon_k)\), with
enlarged polynomial factors in \(\varepsilon_k\).

We use the centered Hessian with a penalty only for the
standard representation in unused degree zero. Write
\(u=(f+g)/\sqrt2\), \(v=(f-g)/\sqrt2\). In active degree,
the proof of
[(VH4)](NEGATIVE_VOTING_HESSIAN.md#the-two-signs-and-the-low-degree-blocks)
gives the antisymmetric lower bound for \(D_A-K_r^0/(r\ell)\)
on the whole space, since the only use of the antisymmetric
mass constraint was to remove the degree-zero term.
For \(D_A+K_r^0/(r\ell)\), removing that term affects only
the standard representation. The standard scalar estimate
in (VH4), omitting the favorable constant term, gives its
gap with the subtraction \(C(\log k)\mathbb E(h)^2\)
in the normalized scalar model of that proof. To convert the
penalty, a standard component associated with a centered
label vector \(a\) has speeds \((a_i-a_j)h\) in the
corresponding facet coordinates. Its common representation
and surface factor is \(k v_k\|a\|^2=\ell\|a\|^2\), while
its mass vector is \(\ell a\,\mathbb E(h)\). Thus the penalty
in the weighted form is at most
\[
 C(\log k)\|Cu\|^2/\ell .
\]
Orthogonality of the representation coordinates gives the same
bound throughout the standard block. All other blocks have
the same bounds as in (VH4). The nonzero unused degrees
retain every term of the interaction and have the bounds from
the dimension lift (NL7). The moment gauge after rotation
correction controls the full norm by a polynomial multiple
of the norm off the common rotations, by the moment-map
inverse above. Consequently for an absolute \(B_0\),
\[
 D_A(f,f)+D_A(g,g)+2K_r^0(f,g)/(r\ell)
 \ge c k^{-B_0}\bigl(D_A(f,f)+D_A(g,g)\bigr)
            - C(\log k)\|Cu\|^2/\ell .
 \tag{PN11}
\]

The exact constant term of the original volume interaction
in units \(r\ell J\) is \(\|m\|^2/(r\ell J)\).
Since \(\|Cu\|^2\le4\|m\|^2/J+O(\varepsilon_k^2)\),
the constant term absorbs the entire \(\|m\|^2/J\) penalty
in (PN11) when \(k\) is large, using \(r\log k=o(1)\).
Thus (PN9)–(PN11) bound the normalized bilinear deficit
below by \(c k^{-B_0}-O(\varepsilon_k)\). First choose
\(\xi=k^{-B_1}\) with a sufficiently large absolute \(B_1\),
then take \(A\) large enough in (PN1) that \(J\) makes
\(\varepsilon_k=o(k^{-B_0})\). The deficit is strictly
positive when \(J>0\), proving (PN1) and its equality case.
