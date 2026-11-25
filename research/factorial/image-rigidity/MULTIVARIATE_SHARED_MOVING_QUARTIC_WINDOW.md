# Uniform moving quartic windows for varying fully shared supports

Let \(a,b\) be distinct ordinary fully shared exponent columns.
Remove equal rows and write \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\),
\(W=\sum_i d_i^2/m_i\), and \(\Delta=\sum_i d_i\).
Use the coefficient phase
\[
 G(q)=\sum_i[(a_i+d_iq)\log(a_i+d_iq)]_{\rm c},\qquad
 F(q)=G(q)-q\log q-(1-q)\log(1-q),
\]
where the chord is through zero and one. Put
\(\rho(q)=q(1-q)G''(q)\) and
\(\delta=\max_{0<q<1}\rho(q)-1\).
For profiles with \(|\delta|\) sufficiently small on a fixed positive
bounded \(W\)-band, let \(q_*\) be the real zero of \(F'''\) near the
maximizer of \(\rho\), and put \(w_*=F'(q_*)\).
The defining neighborhood and uniqueness are established below.
Use the centered parameters at either order \(j=n,n+1\),
\[
 w_j=-\operatorname{Log}(j^\Delta t)-C_0,qquad
 C_0=\sum_i(b_i\log b_i-a_i\log a_i-d_i).
\]

**Theorem.** Fix \(0<W_-\le W_+<\infty\) and \(K<\infty\).
There is a radius \(r_*>0\) and a terminal cutoff such that,
uniformly for all larger \(n\) and all profiles with
\(W_-\le W\le W_+\), \(\sqrt n\,|\delta|\le K\), there is no
consecutive common root \(t\) with \(|w_j-w_*|<r_*\) at either
order \(j=n,n+1\). The logarithms are taken near the positive
real scales in this neighborhood. The degree difference
\(\Delta\) is unrestricted. Neither the radius nor the cutoff
is evaluated. In the bounded quartic window the comparison
uses the imaginary motion of root branches, since the exact
order drift \(w_{n+1}-w_n=-\Delta\log(1+1/n)\) is real.

## Uniform centers and finite integral reduction

Write \(q_m\) for the unique maximizer of \(\rho\). By the
[bounded profile estimates](MULTIVARIATE_VARYING_PROFILE_PHASE.md),
the ratios \(|d_i|/(2m_i)\) are uniformly bounded by some \(R<1\).
The curvature transform and its derivatives are uniformly
bounded in a fixed complex neighborhood of \([0,1]\).
Since \(\rho=q(1-q)G''\), its maximizer stays in a compact
subinterval of \((0,1)\) when \(\max\rho\) is close to one.
In each row, for real \(0<q<1\),
\[
 \rho_i''(q)=-\frac{2a_i b_i d_i^2}{(a_i+d_iq)^3}
      \le -\frac{2(1-R^2)}{(1+R)^3}\frac{d_i^2}{m_i}.
\]
Consequently \(\rho''\) is uniformly negative there. Since
\(F''=(\rho-1)/(q(1-q))\), one has
\(F'''(q_m)=O(\delta)\) and \(F''''(q_m)\le-c<0\) for small
enough \(|\delta|\). Uniform bounds on further derivatives
preserve the latter inequality in a fixed small real
neighborhood of \(q_m\). Strict monotonicity and the intermediate
value theorem give the unique \(q_*\) in this neighborhood,
with \(q_*-q_m=O(|\delta|)\). Define
\[
 a_2=F''(q_*)/2=O(|\delta|),\qquad
 g=F''''(q_*)/24<0,qquad \ell=(-1/(4g))^{1/4}.
\]
The constants \(-g\) and \(\ell\) are uniformly bounded above
and away from zero. For \(\beta=w_*-w\) the real entropy phase
\(\Phi_w(q)=F(q)-wq\) has the exact local form
\[
 \Phi_w(q_*+h)-\Phi_w(q_*)=\beta h+a_2 h^2+g h^4+h^5 E(h),
 \tag{1}
\]
where \(E\) is analytic and uniformly bounded, with its
derivatives on a smaller fixed disk.

In a putative counterexample sequence, extract a retained
and diffuse limit using the same compactness result. Quantities
associated to successive profiles here can vary arbitrarily
slowly with \(n\). In the limit \(\max\rho_*=1\), attained at
a unique \(q_c\), and \(q_*\to q_c\), \(w_*\to w_c=F_*'(q_c)\).
The limiting real function \(F_*\) is strictly concave, since
\(F_*''\le0\) with equality only at \(q_c\).
For \(w=w_c\) the only real maximum is \(q_c\).
For each other fixed real \(w\) it has a unique nondegenerate
real maximum, because the entropy derivative tends to opposite
infinities at the endpoints.

The normalization of \(M_j(t)\) by \(\prod_i(ja_i)!\) gives an
analytic interpolation of its summand in a fixed rectangle
inside \(0<\Re q<1\). Uniform Stirling writes it as
\[
 j^{-1/2} A_j(q)e^{j\Phi_w(q)},\qquad
 A_j(q)=A_0(q)+O(j^{-1}),\qquad
 A_0(q)=\frac{\exp(\tfrac12\sum_i\log((a_i+d_iq)/a_i))}
                  {\sqrt{2\pi q(1-q)}} .                    \tag{2}
\]
Here \(A_0(q_*)\) is positive and uniformly bounded above
and away from zero. The logarithmic amplitude is relative
to the zero endpoint, as required by centering at \(C_0\).
Indeed the terms linear in \(j\) from the row Gamma factors
give \(jG(q)+jq(C_0+\Delta\log j)\), and the latter linear
term cancels against that part of \(t^{jq}\).
These estimates hold for complex \(j\) on a fixed narrow
sector around the positive order ray, including all
\(q\) derivatives on smaller rectangles. In that sector
row arguments are a fixed distance in angle from the
negative real axis and have modulus comparable to
\(|j|m_i\). The relative Stirling logarithmic errors
sum to \(O(W/|j|)\), using \(\sum_i1/m_i\le W\) because
the nonzero gaps are integers. Cauchy's estimate therefore
gives \(\partial_j(A_j-A_0)=O(j^{-2})\) along the positive
real ray.

For \(w\to w_c\), keep a fixed small interval about \(q_*\)
with endpoints first placed at neighboring half-grid points
at order \(j\). At the limiting cusp the real phase has
a strict gap outside any such interval. The ordinary entropy
upper bound for the binomial coefficients and the bounded
Stirling log amplitude on \([0,1]\) give an exponential loss
for the excluded summands, uniformly for the approaching
profiles. Apply the
[cotangent summation argument](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour)
on the retained interval. The endpoint walls have cotangent
modulus at most one and keep a strict phase gap at a fixed
small height. The horizontal sides have constant cotangent
terms up to \(O(e^{-2\pi j\nu})\) at heights \(\pm\nu\).
The phase increase there is \(O(\nu^2)+|\Im w|\nu\)
by real analyticity. Choose \(\nu\) small, then pass sufficiently
far along the sequence to absorb this increase. The two
constant-kernel integrals deform to the real interval and
add. Moving the half-grid endpoints to the fixed coordinates
\(q_*\pm r_0\) costs another exponential remainder for some
small fixed \(r_0>0\). Thus, after the scalar normalization,
the exact value is
\[
 j^{1/2}\int_{q_*-r_0}^{q_*+r_0} A_j(q)e^{j\Phi_w(q)}\,dq
       +O\bigl(e^{j(\Re\Phi_w(q_*)-c)}\bigr)                 \tag{3}
\]
for uniform \(c>0\) on a small neighborhood of \(w_c\).
For a fixed real \(w\ne w_c\), the same argument on an interval
around its limiting real maximum gives the usual uniform
one-saddle Gaussian expansion and local eventual zero-freeness,
since that stationary point is nondegenerate.

## Bounded quartic fields and the adjacent imaginary drift

For a given profile and positive real order parameter set
\[
 \lambda=j^{-1/4},\qquad q=q_*+\ell\lambda y,qquad
 T=-2a_2\ell^2/\lambda^2,qquad H=\ell\beta/\lambda^3 .
\]
The variables \(T,H\) can be treated as independent in the
local interpolation, by inserting the quadratic and linear
terms they prescribe into (1) and keeping the remainder
and amplitude fixed. After dividing the integral in (3)
by its centered exponential and by
\(\ell\lambda A_0(q_*)\), the resulting function
\(J(T,H,\lambda)\) tends to
\[
 Z(T,H)=\int_{\mathbb R}e^{-y^4/4-Ty^2/2+Hy}\,dy .
 \tag{4}
\]
The convergence holds uniformly with the \(T\)- and
\(H\)-derivatives for bounded real \(T\) and bounded complex
\(H\), and \(\partial_\lambda J\) is uniformly bounded
there. In fact the remaining exponent is
\(\lambda\ell^5 y^5 E(\ell\lambda y)\), and by reducing
\(r_0\) the real part of the exponent is bounded by
\(-c_1y^4+C y^2+C|y|\) over the entire scaled interval.
Each differentiated factor is bounded by a polynomial
in \(|y|\) under this envelope. For the amplitude's
order derivative, (2) gives
\(\partial_\lambda A_{\lambda^{-4}}=O(\lambda^3)\) at
fixed \(q\). The moving endpoints at \(y=\pm r_0/(\ell\lambda)\)
are exponentially small. These estimates prove the
differentiated convergence as well as
\(J-Z=O(\lambda)\) and permit differentiation on the
real order ray. At actual integer orders, the normalized
error in (3) is exponentially small uniformly on bounded
quartic fields.

For real \(T\), the zeros of \(Z(T,\cdot)\) are simple and
purely imaginary, with positive imaginary parts denoted by
\(y_k(T)\) in increasing order, and their negatives.
The
[quartic root-speed inequality](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#a-quartic-root-speed-inequality)
and its
[continuation to both signs of the temperature](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#the-bounded-moving-critical-field)
give
\[
                       3y_k(T)-2T y_k'(T)>0 .              \tag{5}
\]
Only finitely many model branches enter a given bounded
range for compact \(T\). These branches are uniformly
separated, and their derivatives in (5) have a uniform
strict margin for each fixed finite set of indices.
Uniform convergence and simplicity therefore give
corresponding complex zero branches \(H_k(T,\lambda)\)
for the local interpolation. For the branch at positive
imaginary part, \(H_k\to i y_k(T)\) and
\(\partial_T H_k\to i y_k'(T)\), while
\(\partial_\lambda H_k\) is bounded. The negative-imaginary
branches have conjugate limits. These claims follow
from the implicit function theorem and the uniform
bounds just proved. Using a slightly larger field
disk, Rouché's theorem also puts the actual integer-order
zeros exponentially close to these branches.

With the profile fixed across adjacent orders, use
\(T(\lambda)=-2a_2\ell^2/\lambda^2\). Along a branch its
physical field is
\[
 \beta_k(\lambda)=\ell^{-1}\lambda^3 H_k(T(\lambda),\lambda),
 \qquad
 \beta_k'(\lambda)=\ell^{-1}\lambda^2
         (3H_k-2T\partial_T H_k+\lambda\partial_\lambda H_k).
\]
By (5) the imaginary part of this derivative has strict
constant sign and modulus at least \(c_2\lambda^2\)
on each selected branch, uniformly for bounded temperatures.
Thus the two adjacent physical imaginary parts on the
same branch differ by at least a positive multiple of
\(n^{-7/4}\). The exponentially small root errors cannot
remove this drift. At bounded quartic fields, roots on
different branches cannot have equal physical imaginary
parts either. Indeed \(\Im H\) at a common original \(t\)
can only change by the factor \((1+1/n)^{3/4}\), while the
temperatures differ by \(O(n^{-1})\) and distinct limiting
imaginary coordinates have a fixed separation.
This excludes a common root whenever both scaled fields
are bounded, for any real order drift in the physical
parameter.

## Escaping fields at bounded scaled temperature

To prove the asserted fixed-radius statement, it suffices
to exclude sequences of common roots with one of
\(w_n,w_{n+1}\) tending to \(w_*\) in difference. Otherwise
one could take the radius to zero and the terminal to
infinity along a counterexample sequence. Extract a
retained and diffuse limit as above. The other centered
scale stays in a compact set by the
[uniform root annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md)
and the [centered phase estimates](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#centering-compact-scales-and-uniform-local-expansions).
Its imaginary part is unchanged. A limit away from
\(w_c\) would be a noncusp real scale and is excluded by
the one-saddle estimate following (3). Thus both fields
\(\beta=w_*-w_n\) and \(\beta'=w_*-w_{n+1}\) tend to zero.

Consider one order \(j\) whose field satisfies
\(j^{3/4}|\beta|\to\infty\), denoting that order's
field by \(\beta\) for the moment. Let
\[
 \epsilon=|\beta|^{1/3},\quad
 \ell_0=(-1/(4g))^{1/3},\quad h=\ell_0\epsilon z,
 \qquad M=j\ell_0\epsilon^4\longrightarrow\infty .
\]
After dividing the centered phase (1) by
\(\ell_0\epsilon^4\), its scaled form is
\[
 \Psi(z)=e^{i\arg\beta}z-z^4/4
                   +(a_2\ell_0/\epsilon^2)z^2+O(\epsilon|z|^5).
 \tag{6}
\]
The coefficient of \(z^2\) tends to zero, since
\(\sqrt j\,|a_2|\) is bounded and
\(\sqrt j\,\epsilon^2\to\infty\). The convergence on
bounded disks includes the local phase derivatives.
At bounded heights, the large-real-part tail on
the retained physical disk has real part at most
\(-c_3|\Re z|^4\), by the small choice of \(r_0\) in
(3). Indeed \(\epsilon|z|^5\) is at most a small
multiple of \(|z|^4\) there.

Use the [horizontal quartic contours](MULTIVARIATE_SHARED_REAL_CUSP_ESCAPING_FIELDS.md#one-saddle-and-two-saddle-contours)
for the limiting phase in (6). When the limiting
argument is separated from \(\pi/2,-\pi/2\), there is
a simple stationary point on whose horizontal line
the real phase has a unique strict maximum.
To recall the check, for \(0\le\theta<\pi/2\) this
is \(z=e^{i\theta/3}=u+iv\), and the real phase on
its line is a constant plus
\(-x^4/4+3v^2x^2/2+x\cos\theta\).
Positive \(x\) wins over its reflection, and the
derivative has exactly one positive zero. Its
curvature there is negative. Reflection and conjugation
cover the other angles away from the transitions.
The corresponding actual horizontal contour for (6)
has the same strict loss and the quartic tail.
Its Gaussian expansion has nonzero amplitude.
The fixed physical remainder in (3) is negligible
relative to that packet since \(M=o(j)\).
Thus every escaping zero must have field angle
tending to \(\pi/2\) or \(-\pi/2\).

If one of two common-root fields is bounded in the
quartic scale, the imaginary part of the other is
also physically \(O(n^{-3/4})\). An escaping other
field would have an asymptotically real angle, which
was just excluded. When both fields are bounded, use
the preceding bounded-field comparison. It remains
to exclude both fields escaping. They share the same
imaginary part. By conjugation take that imaginary
part positive. Then both arguments tend to \(\pi/2\),
and
\[
            \beta'/\beta\longrightarrow1,\qquad
            \beta'-\beta=\Delta\log(1+1/n)\in\mathbb R .
 \tag{7}
\]

For angles sufficiently close to \(\pi/2\), use the
two half-contours through the actual saddles near
\(\pm\sqrt3/2+i/2\) in (6). The horizontal line through
each limiting saddle has a unique strict maximum on its
half-line. The two real stationary heights are \(-3/8\)
and the midpoint height is \(-33/64\). The midpoint
connectors therefore have a uniform strict scaled loss
after the small perturbations in (6). The connectors
at the fixed physical endpoints retain their physical
loss from (3). In particular the two actual packets
give uniform separate Gaussian expansions with relative
error \(O(M^{-1/2})\) and nonzero bounded leading amplitudes
after extracting the common length factor
\(\epsilon M^{-1/2}\) and their respective stationary
exponentials. Here and below take \(\epsilon,M\) from the
lower order when both orders occur in one expression.
The analogous estimates at the upper order use quantities
with asymptotically the same sizes, by (7). After the
common centered exponential is removed, the fixed
physical remainders are \(O(e^{-c n})\) and the individual
packet moduli are at least
\(c_4\epsilon M^{-1/2}e^{-C M}\).

Write \(q_R(\beta),q_L(\beta)\) for the two physical
saddles, and let \(D(\beta)\) be their right-minus-left
stationary phase difference. Modulus cancellation at
a root and the separate expansions imply
\[
            \Re D(\beta)=O(1/n),qquad
            \Re D(\beta')=O(1/n).
 \tag{8}
\]
The saddle branches are analytic back to each other on
the segment joining the fields in (7), using the bounded
scaled chart at the lower order. Stationarity and the
linear field in (1) give
\[
 D_\beta=q_R-q_L,qquad
 \Re(q_R-q_L)=\ell_0\epsilon(\sqrt3+o(1))
 \tag{9}
\]
uniformly on that segment. Integrate (9) on the real
displacement in (7) and use (8). This gives
\(|\Delta|\epsilon=O(1)\).

Retain exactly the two lower-order contours at the upper
order. After removing the centered exponential at each
order, their integrand ratio is \((A_{n+1}/A_n)e^{B_n(q)}\),
where the linear dependence on \(w\) gives
\[
 B_n(q)=\Phi_{w_n}(q)-\Phi_{w_n}(q_*)
               +\Delta(n+1)\log(1+1/n)(q-q_*).
 \tag{10}
\]
The amplitude ratio is \(1+O(n^{-1})\) on the fixed
physical disk. Near a lower saddle \(q_j\), at its scaled
coordinate \(z_j\), (10) gives
\[
 |B_n(q)-B_n(q_j)|
          \le C\bigl(\epsilon^4 |z-z_j|^2
                            +|\Delta|\epsilon |z-z_j|\bigr),
]
by exact lower-order stationarity. Over a whole retained
contour its bound is
\(C[\epsilon^4(1+|z|^4)+|\Delta|\epsilon(1+|z|)]\).
Consequently after extracting \(e^{B_n(q_j)}\), replacing
the remaining exponential by one gives a relative error
\[
        O\left(\frac{\epsilon^4}{M}
          +\frac{|\Delta|\epsilon}{\sqrt M}
          +\frac1n\right)
 \tag{11}
\]
for that packet. The linear and quadratic local factors
give their indicated Gaussian moments. The bounded
scaled pieces away from the saddle have strict phase
loss and absorb polynomial powers of \(M\). On the
remaining horizontal tail, \(|\Delta|\epsilon=O(1)\)
makes the additional exponential grow at most like
\(\exp(C(1+|z|)+C\epsilon^4(1+|z|^4))\), absorbed by
the quartic majorant. The physical connectors with
their fixed gap absorb the multiplier growth there,
which is at most \(e^{O(|\Delta|+1)}=e^{o(n)}\).
This proves (11) using the same contour for the two
orders. All comparisons refer to the centered
normalizations so that the multipliers at the saddles
stay bounded above and away from zero.

If \(\Delta\ne0\), the real difference of the exponents
in (10) at the two saddles is
\[
 \Re D(\beta)+\Delta(n+1)\log(1+1/n)\Re(q_R-q_L).
\]
Its magnitude is at least a positive multiple of
\(|\Delta|\epsilon\) by (8) and (9), since
\(1/n=o(\epsilon^4)\). The two exponential multipliers
therefore differ in modulus by at least that order.
If \(\Delta=0\), their exponent difference is
\(D(\beta)\), whose imaginary part is asymptotic to
\((3\sqrt3/4)\ell_0\epsilon^4\) in the model (6).
Thus their difference in modulus as complex numbers
has order at least \(\epsilon^4\). In both cases (11)
is smaller by a factor tending to zero. At a common
root the two packets would cancel at both orders,
forcing their two adjacent ratios to agree up to
exponentially small errors. The physical remainders
stay negligible at these precisions since \(M=o(n)\),
\(n\epsilon^4\to\infty\), and the packets have the
lower bounds above. This contradiction excludes
both escaping fields and proves the theorem.
The original coefficients are integers, so rational
Bézout gives the characteristic-zero assertion. \(\square\)
