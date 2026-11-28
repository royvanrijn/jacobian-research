# A radial dual for variable angular slopes on finite shells

Take any finite partition of Gaussian radii into consecutive
positive-mass shells indexed by \(i\). Consider continuous,
even angular kernels \(F_{\varepsilon,ij}\) forming a one-sided
perturbation of the regular three-cell simplex agreement
\(F_T(t)=(1-t/a)_+\), \(a=2\pi/3\), \(0\le t\le\pi\).
Require the range, triangle, balance and four-point constraints.
Assume the difference quotients are uniformly bounded and
their derivatives \(f_{ij}(t)\) exist away from \(a\), with
\[
  f_{ii}(0)=0,\qquad
  f_{ii}(t)=A_i t/a+o(t),\qquad
  f_{ij}(t)=f_{ij}(0)+o(1)\quad(t\downarrow0).
\]
The slopes \(A_i\) may differ and can have either sign.
Write \(b_{ij}=-f_{ij}(0)\), and let \(r_j\) be the radius
between shells \(j,j+1\). Define the adjacent radial crossing
cost by \(C=\sum_j b_{j,j+1}r_j e^{-r_j^2/2}\). Put
\(A(r)=A_i\) on shell \(i\). For nonnegative slopes,
\[
  G:=\int_0^\infty \frac{A(r)}{a r}\,r e^{-r^2/2}\,dr
                    \le \frac{997}{1000} C.                     \tag{VR1}
\]
For signed slopes the constant in (VR1) becomes \(999/1000\),
as proved below. In both cases the first variation of the relaxed Mehler
objective as \(\rho\uparrow1\) has expansion
\[
   \left.\frac{d}{d\varepsilon}\mathcal S_\rho(F_\varepsilon)
           \right|_{0+}
       =\frac2{\sqrt\pi}(G-C)\sqrt{1-\rho}
                  +o(\sqrt{1-\rho}).
\]
If any slope or adjacent jump is nonzero, the leading
coefficient is strictly negative. The result concerns the
first variation for the stated finite-shell perturbations.
[The common-slope case and the saturated clique](GAUSSIAN_AGREEMENT_TRIANGLE_STRENGTHENINGS.md#a-common-angular-slope-across-finitely-many-radial-shells)

## Row and four-point bounds for varying slopes

Use the radial probability coordinate \(u=1-e^{-r^2/2}\),
write \(r(u)=\sqrt{-2\log(1-u)}\) and identify \(A,b\)
with their step functions in the probability coordinate.
Let \(\ell=\sum_j b_{j,j+1}\delta_{u_j}\) at the shell
boundaries. Active triangles at zero make \(b\) a
pseudometric. Consequently \(b(u,v)\) is at most the
\(\ell\)-length between \(u\) and \(v\).
At \(0<t<a\), subdivide active angular triangles with
intermediate points in the left endpoint shell, then with
intermediate points in the right endpoint shell. The two
inequalities give
\[
 f_{ij}(t)\ge \max(A_i,A_j)t/a-b_{ij}.
\]
For \(t>a\) the derivative is nonnegative by the range
constraint. The exact row mean therefore gives
\[
 \int_0^1\max(A(u),A(v))\,dv\le 2 B(u),\qquad
 B(u):=\int_0^1 b(u,v)\,dv.                                     \tag{VR2}
\]
All angular evaluations of the constraints, including at
coincident angles, are justified by continuity of each
perturbed kernel. Differentiation of balance uses the uniform
bound on the difference quotients.

Set \(m(c)=\min(c,1-c)\), \(d(c)=m(c)|c-1/2|\),
\(M=\int m\,d\ell\), \(D=\int d\,d\ell\) and
\(\bar A=\int_0^1 A\). The median row and path bound give
\(\bar A\le2M\).
Let \(T_0\) be the mean of \(\max(A(u),A(v))\) for two
independent uniform probabilities. Let \(T_\times\) be its
mean conditional on opposite median halves. Write
\(\bar b,b_\times\) for the corresponding means of \(b\).
Subtract the angular baseline above on \([0,a)\) and zero
on \((a,\pi]\). The correction is nonnegative and its total
integral over radial pairs and \(0<t<\pi\) is
\(a(\bar b-T_0/2)\) by balance. The alternating median
four-point test used in (TS7) has baseline expectation
\(3T_\times-4b_\times\). Its angular and radial density
bounds give a correction at most \(16(\bar b-T_0/2)\).
Thus
\[
 8T_0-3T_\times\le 16\bar b-4b_\times\le8M+32D .
\]
For the last inequality use the path-length estimates (TS6).
Since \(T_0=\bar A+\tfrac12\mathbb E|A(u)-A(v)|\) and
\(T_\times=\bar A+\tfrac12\mathbb E_\times|A(u)-A(v)|\),
with \(\mathbb E_\times|A(u)-A(v)|\le2\mathbb E|A(u)-A(v)|\),
one obtains
\[
            5\bar A+\mathbb E|A(u)-A(v)|\le8M+32D.              \tag{VR3}
\]

## A weighted allocation of the row bounds

Write \(\omega(u)=1/(a r(u))\), choose \(\alpha=29/200\)
and \(\beta=3/40\), and define
\[
 f(u)=(\omega(u)-\alpha)_+,\qquad Z=\int_0^1 f,\qquad
 J(c)=\frac1{1-c}\int_c^1 f,\qquad \lambda=-J'.
\]
The weight \(f\) is integrable and decreasing. Thus \(\lambda\)
is nonnegative, \((1-u)\lambda(u)+J(u)=f(u)\) almost
everywhere, and \(\int_c^1\lambda=J(c)\) since \(f\) vanishes
on a neighborhood of one. In particular \(\int_0^1\lambda=Z\).

In (VR2) weight the row at \(u\) by \(\lambda(u)\) and
allocate the maximum for \((u,v)\) to the endpoint of smaller
radial probability. Integration gives the coefficient
\((1-u)\lambda(u)+\int_u^1\lambda=f(u)\) on \(A(u)\).
For a radial cut at \(c\), the weighted path cost is
\((1-c)(Z-J(c))+cJ(c)\), by counting independently placed
endpoints on opposite sides of the cut. Hence
\[
 \int_0^1 f(u)A(u)\,du
         \le2\int[(1-c)(Z-J(c))+cJ(c)]\,d\ell(c).
\]
Since \(A\ge0\), the pointwise inequality \(\omega\le f+\alpha\)
gives \(G\le\int fA+\alpha\bar A\). Price
\((\alpha-\beta)\bar A\) by the median row and \(\beta\bar A\)
by (VR3), dropping its nonnegative difference term. It follows
that \(G\le\int W(c)\,d\ell(c)\), with
\[
 W(c)=2[(1-c)(Z-J(c))+cJ(c)]+2\alpha m(c)
                   +\beta[-2m(c)/5+32d(c)/5].                  \tag{VR4}
\]
The radial density satisfies \(w(c)=(1-c)r(c)\) and
\(C=\int w(c)\,d\ell(c)\).

## Exact one-dimensional bound and the Gaussian limit

The [rational interval checker](check_agreement_four_point_radial_dual.py)
proves \(W(c)<(997/1000)w(c)\) for \(0<c<1\). It checks
the radius interval \([1/10,33/10]\) by subdivision with
rational Taylor enclosures of \(e^{-r^2/2}\) and
\(\int_0^r e^{-s^2/2}\,ds\). It uses Machin's formula with
alternating arctangent bounds for \(\pi\). In radius
coordinates the primitive of \(f\) up to its cutoff
\(1/(a\alpha)\) is
\[
  \frac1a\int_0^r e^{-s^2/2}\,ds-\alpha(1-e^{-r^2/2}).
\]
Beyond the cutoff the primitive is the constant \(Z\).
The checker uses monotonicity of \(J\) and of this truncated
primitive. For the left half, with \(K(c)=\int_0^c f\) and
\(s=1-c\), the angular term inside the factor two in (VR4)
is \([c^2 Z+(1-2c)K(c)]/s\). This also bounds the small
radius tail using \(c\le r^2/2\), \(K(c)\le r/a\) and
\((1-2c)/s^2\le1\). The remaining terms in (VR4) are at
most \(m(c)/2\). For \(r\ge33/10\), \(J=0\) and
\(W(c)/w(c)\le(2Z+1/2)/r\). Both tail bounds are strictly
below \(997/1000\). Integrating the certified bound relative
to \(\ell\) proves (VR1).

For the objective expansion, given \(X\) with \(|X|=r>0\),
the angular gap \(\vartheta\) to
\(Y=\rho X+\sqrt{1-\rho^2}Z'\) satisfies
\(\vartheta/\sqrt{1-\rho}\to\sqrt2|Z'_2|/r\) almost surely
after choosing the radial axis along \(X\).
For \(\rho\ge1/2\) this quotient is bounded by
\(2\pi\sqrt2|Z'|/r\), using the Euclidean noise displacement
and the angle to its endpoint. The bound is integrable in
the planar Gaussian radial measure. Dominated convergence
gives the coefficient \(2G/\sqrt\pi\) from same-shell
angular slopes. Discarding radial crossing events from
that angular term changes it by \(o(\sqrt{1-\rho})\),
by the angular tail estimate and crossing estimate in (TS4).
The same estimates control the uniform small-angle
remainders over a finite number of shells. Each adjacent
boundary crossing contributes its radial Gaussian
perimeter times \(2b_{j,j+1}\sqrt{1-\rho}/\sqrt\pi\)
to the radial loss. Nonadjacent crossings require
positive minimum Euclidean displacement and have
exponentially small probability. This proves the
expansion and the strict sign under the stated condition.
[Angular and radial crossing estimates, TS4](GAUSSIAN_AGREEMENT_TRIANGLE_STRENGTHENINGS.md#consequence-for-the-near-one-objective)

## Signed angular slopes

Under the same hypotheses with arbitrary real \(A_i\), one has
\(G\le(999/1000)C\). The near-one expansion above is unchanged.
Its leading coefficient is strictly negative if any adjacent
jump or angular slope is nonzero.

The row and clique inequalities preceding (VR3) hold for
signed \(A\), since they use the pointwise maximum of two
slopes. The allocation in (VR4) would leave an excess
coefficient on \(A\) in the far radial tail. We remove that
coefficient by reallocating the maxima in the median row and
the clique bound. Keep the constants \(\alpha=29/200\),
\(\beta=3/40\), and set \(\mu=\alpha-\beta=7/100\).
Let \(u_\alpha\) solve \(\omega(u_\alpha)=\alpha\), and
choose \(u_0<u_\alpha\) by
\(\omega(u_0)=\alpha+1/1000=73/500\). Write \(E=[u_0,1]\)
and \(\tau=1-u_0\).

Use the median row contribution with weight \(\mu\),
averaged over the uniform law \(\nu\) on \([12/25,13/25]\).
In its maximum for a row sampled from \(\nu\) and an
independent uniform endpoint \(v\), allocate to \(v\) outside
\(E\). For \(v\in E\), allocate to \(v\) with fraction
\(\min(1,\omega(v)/\alpha)\) and to the row endpoint
otherwise. These fractions give a linear lower bound on the
maximum for either sign of each slope. The resulting
coefficient in \(E\) is
\(g_M(v)=\mu\min(1,\omega(v)/\alpha)\).
Off \(E\) it is \(\mu\) plus a bump of mass
\(\mu H/\alpha\) distributed according to \(\nu\), where
\[
                       H=\int_{u_\alpha}^1(\alpha-\omega).
\]

For the clique use its full left-hand side
\(8T_0-3T_\times=T_\times+2T_{--}+2T_{++}\) with weight
\(\beta/5\). Here \(T_{--},T_{++}\) are the means of the
maximum on the two respective same-half products. Allocate
equally between endpoints outside \(E\). For two endpoints
in \(E\) allocate to the smaller one. On a pair with just
one endpoint \(v\in E\), allocate to it with fraction
\[
 q(v)=
 \frac{\omega(v)-g_M(v)-(16\beta/5)(1-v)}
      {\,2\beta-(16\beta/5)\tau\,}.
                                                               \tag{VR5}
\]
The contribution on \(E\) from \(E\times E\) is
\((16\beta/5)(1-v)\), by integrating the larger endpoint in
the \(2T_{++}\) term. The available coefficient from partners
outside \(E\) is \(2\beta-(16\beta/5)\tau\), consisting of
\(2\beta/5\) from \(T_\times\) and
\((16\beta/5)(1/2-\tau)\) from \(2T_{++}\). Thus the
combined median and clique coefficient is exactly \(\omega\)
on \(E\). On the two halves off \(E\), the clique
coefficient differs from \(\beta\) by constants
\[
 \delta_-=(4\beta/5)\int_E(1/2-q(v))\,dv,\qquad
 \delta_+=4\delta_- .
                                                               \tag{VR6}
\]
These formulas follow by assigning the unused fraction
\(1-q(v)\) to the other endpoint. The densities on the
conditional opposite-half and same-half products are both
four, with two orders for the same-half contribution.

We give explicit bounds for these allocations. Using
\(3.14<\pi<22/7\) gives \(r(u_0)>163/50\) and hence
\(\tau<1/200\). For the latter, the first thirteen
terms of the exponential series at \((163/50)^2/2\)
already exceed \(200\). The numerator in (VR5) is
nonnegative since \(\omega-g_M\ge(\beta/\alpha)\omega\)
and \((16\alpha/5)(1-v)\le\omega(v)\) on \(E\).
The last inequality holds at \(u_0\) by \(\tau<1/200\)
and extends to the tail because \(r e^{-r^2/2}\)
decreases for \(r\ge r(u_0)>1\).
The numerator is at most \(73/500\), while the
denominator exceeds \(2\beta-(16\beta/5)/200>73/500\).
This proves \(0\le q\le1\).

Integrating the two available clique contributions on
\(E\), with fraction \(1/2\) on the outside partners,
gives \(\beta\tau\). The allocated clique coefficient
\(\omega-g_M\) has integral
\[
 \beta\tau+\int_{u_0}^{u_\alpha}(\omega-\alpha)
                                      -(\beta/\alpha)H .
\]
The positive correction is less than \(10^{-3}\tau<5\cdot
10^{-6}\). On \(r\ge7/2\) one has \(\omega<137/1000\)
and Gaussian radial mass \(e^{-49/8}>1/500\). The latter
follows from \(e<11/4\) and
\(e^{1/8}\le8/7\). Therefore \(H>16\cdot10^{-6}\) and
\(\beta/\alpha>1/2\), so the allocated mass is below
\(\beta\tau\). This gives \(\delta_-\ge0\) and
\(\delta_+\ge0\), with \(\delta_+<(8\beta/5)/200=0.0006\).

Let \(\tilde f\) be the coefficient still needed after
these two allocations. It is zero on \(E\) and, on the
two halves off \(E\), respectively equals
\[
 \tilde f(u)=\omega(u)-\alpha-\delta_\pm
                  -\frac{\mu H}{\alpha}\frac{d\nu}{du}.
\]
Then \(0\le\tilde f\le f\). Indeed outside the median band
and \(E\) the uncorrected coefficient is at least \(0.001\),
and in the median band the correction is at most
\(0.0006+\mu\tau/(1/25)<0.01\), much less than the
uncorrected coefficient. Also \(\|\tilde f-f\|_\infty<0.01\),
including \(E\) where \(f\le0.001\).
Define \(\tilde Z,\tilde J,\tilde\lambda\) from \(\tilde f\)
as \(Z,J,\lambda\) were defined from \(f\). Above \(13/25\)
the function \(\tilde f\) is decreasing up to its downward
jump to zero, hence \(\tilde f\ge\tilde J\).
At \(u\le13/25\), the integral over \(v\in[4/5,1]\)
in \(f(u)-J(u)\) gives
\[
 f(u)-J(u)\ge
 \tfrac15[\omega(13/25)-\omega(4/5)]>0.02 .
\]
For the last bound use \(r(13/25)<5/4\) and
\(17/10<r(4/5)<2\), from the exponential series and
\(e<11/4\), then \(a<44/21\).
Since the changes in both \(f\) and \(J\) are less than
\(0.01\) in absolute value, this proves
\(\tilde f\ge\tilde J\) everywhere almost surely. Thus
\(\tilde\lambda=-\tilde J'\) is nonnegative and
allocating its row maxima to the smaller radial
probability gives the exact residual coefficient \(\tilde f\).
The three allocations now sum to \(\omega\) even when
the slopes are signed.

The weighted-row cut cost for \(\tilde f\) is no larger
than for \(f\) at any \(c\). For \(c\ge1/2\) this follows
directly from
\((1-c)\tilde Z+(2c-1)\tilde J(c)\) and the pointwise
inequality \(\tilde f\le f\). For \(c<1/2\) use the form
\([c^2\tilde Z+(1-2c)\int_0^c\tilde f]/(1-c)\), with
both numerator coefficients nonnegative. Reallocating the
clique maxima does not change their cost from (VR3).
Spreading the median row changes its cost only for
\(c\in[12/25,13/25]\). For \(\eta=1/50\) and
\(h=|c-1/2|\le\eta\), the increase is
\(2\mu h(\eta-h)/\eta\le\mu\eta/2=0.0007\).
On this band \(w(c)>12/25\), since \(r(c)>1\).
The original \(997/1000\) certificate now gives a
total cost less than \((999/1000)w(c)\) at every cut.
This proves the signed bound with the same checker
and the stated margin.

The leading objective expansion uses only the signed
slopes and their boundedness on finitely many shells.
If \(C>0\) the signed bound gives strict negativity.
If \(C=0\), then \(b=0\) by the adjacent path bound
and the row constraint forces \(A_i\le0\) for all \(i\).
A nonzero slope then makes \(G<0\), as required.
