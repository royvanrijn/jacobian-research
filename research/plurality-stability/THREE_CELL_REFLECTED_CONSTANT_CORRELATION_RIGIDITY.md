# Correlation derivative asymptotics on a reflected maximizing graph

Let \(A\) be a balanced planar three-cell global maximizer at \(0<r<1\)
with a central reflection interchanging labels 1 and 2. Use the
[monotone graph coordinates](THREE_CELL_REFLECTION_REDUCTION.md)
in which the right \(3/1\) interface is \((x,h(x))\), \(x>0\).
Write
\[
 F=T_r(1_{A_3}-1_{A_1}),\qquad
 g(x)=\left.\partial_\rho T_\rho(1_{A_3}-1_{A_1})
                   (x,h(x))\right|_{\rho=r}.
\]
The derivative holds the partition and the observation fixed.
If \(g\) is bounded on either side on a tail, then the graph
has an affine asymptote of slope \(m>0\) and \(g\) converges
at a Gaussian rate to the affine value \(g_\infty\) in (RS5).
Write \(t\) for the longitudinal coordinate along this
asymptote. Suppose the pressure \(\beta\ne0\). In the cases
\(m\le1\) or \(\beta>0\), the convergence cannot satisfy
\[
       |g-g_\infty|\le C\exp\!\left[-\frac{r^2t^2}{2\sigma^2}
                                                 +C' t\right]
               \qquad\text{on a tail}                     \tag{CR0}
\]
for any constants \(C,C'\), where \(\sigma^2=1-r^2\).

In the special case of constant \(g\), this gives another
proof of the rigidity in the existing reflected contact
analysis. The differential equation below uses tail data
and extends the restriction to (CR0).
[Earlier proof by continuation of the score level](THREE_CELL_FIRST_CONTACT_RIGIDITY_TESTS.md#reflected-ties-have-at-most-second-order-contact)

In particular, suppose \(A\) ties the simplex \(T\) at \(r\) and
distinct correlations of simplex optimality accumulate at \(r\).
For the fixed-partition deficit \(D_\rho=S_\rho(T)-S_\rho(A)\),
the conditions \(D'(r)=D''(r)=0\) force the simplex in this
reflected class. By (G10), they make \(g\) constant via
every compact zero-flux speed on the right interface.
[The mixed contact inequality, G10](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

## The two possibilities for the derivative at infinity

Use the smooth discretely homogeneous tail profile \(H\) and its
pressure correction (TP1). In particular \(m(x)=H'(x)\) is
periodic in \(\log x\), \(h=H+O(1)\) and \(h'=H'+O(1/x)\).
At the conditional source \(S=rx+\sigma Z\), (TP1) gives
\[
 \frac{rh(x)-h(S)}{\sigma}
       =z_\beta\sqrt{1+m(x)^2}-m(x)Z
                         +O\bigl((1+|Z|+Z^2)/x\bigr)
\]
on \(S\ge a x\) for any fixed \(0<a<r\), where
\(z_\beta=\Phi^{-1}((1+\beta)/2)\).
The omitted probability is Gaussian-small. Since \(\phi'\)
is bounded, the expectation \(A_x\) in (RS6) is
\(\phi(z_\beta)/\sqrt{1+m(x)^2}+O(1/x)\).
Consequently (RS6) gives
\[
 g(x)=
   \frac{2\phi(z_\beta)}{\sigma\sqrt{1+H'(x)^2}}
                   \bigl(H(x)-xH'(x)\bigr)+O(1).          \tag{CR9}
\]
If \(H\) is nonlinear, \(H(x)/x\) is a nonconstant smooth
periodic function of \(\log x\). Its logarithmic derivative
has both signs, so (CR9) yields linearly growing positive
and negative values of \(g\) on two subsequences. A tail
bound on either side therefore forces linear \(H\).
[The pressure-dependent phase correction, TP1](THREE_CELL_REFLECTED_TAIL_PHASE_REGULARITY.md)
[The correlation derivative, RS6](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#a-constant-correlation-derivative-constrains-the-tail)

With \(H(x)=mx\), (TP1) gives \(h(x)=mx+c+O(1/x)\) with
\(c=\sigma z_\beta\sqrt{1+m^2}/(r-1)\) and
\(h'(x)=m+O(1/x)\). In the rotated coordinates used below,
the right interface is a bounded perturbation of its
pressure line with offset \(b\), and \(u(t)\to0\).
Equations (CR2)–(CR3) give a binary pressure equation
with an additive Gaussian-small error throughout a fixed-width
tube about that line, including when \(m>1\). The proof
of (TF1)–(TF6) uses only a Gaussian bound on this extra
error and gives \(|u(t)|\le C e^{-c_0t^2}\). In (TF1)
the error from (CR3) is added to the left Gaussian error,
and the improved residual and iteration have the same form.

Applying the correlation derivative to (CR2)–(CR3) preserves
some Gaussian rate. Its kernel is the conditional density
times a polynomial in \(t,z,s,v\). Integrating that density over
the wall perturbation of width \(O(e^{-c_0 s^2})\), the left
half-plane error, and the opposite wedge gives a Gaussian
bound in \(t\). The difference between evaluating the affine
derivative at \(b+u(t)\) and at \(b\) is Gaussian-small as well.
Thus \(g\) converges at a Gaussian rate to the derivative of
the affine score at \(b\).
[Gaussian flattening and its iteration, TF1–TF6](THREE_CELL_TAIL_FLATNESS.md)

## The pressure line and the second-order equation

The stationary score \(F\) is constant on the graph, say \(\beta\).
In the proof of the restriction (CR0), suppose \(\beta\ne0\)
and the phase is linear. For the constant-\(g\) corollary,
the estimates (RS5)–(RS7) give the required slope case
when \(\beta\ne0\). Namely \(m<1/\sqrt3\) if \(\beta<0\),
and \(m>1/\sqrt3\) if \(\beta>0\). The slope comparison
uses the sign of \(xh'(x)-h(x)\) from the resolvent (RS7)
and balance. When \(g\) is constant and \(\beta=0\),
(RS7) already gives the simplex.
[Constant derivative tail estimates, RS5–RS7](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#a-constant-correlation-derivative-constrains-the-tail)

Take orthonormal coordinates \((t,z)\) along and normal to the
right asymptote, oriented so that the third label is above it.
Write the right tail as \(z=b+u(t)\) with \(b=c/\sqrt{1+m^2}\),
\(u,u'\to0\), and bounded derivatives. In these coordinates
\(t=(x+my)/\sqrt{1+m^2}\) on the original plane. The derivative
of \(t\) along the right interface is positive and bounded away
from zero, so the tail estimates pass to the new graph.
Set \(\sigma^2=1-r^2\), \(\psi(v)=2\Phi(v/\sigma)-1\), and
\(G(z)=\psi(rz-b)\). This is the score of the binary halfspace
contrast with dividing line \(z=b\). Then \(G(b)=\beta\),
and \(b\) has the opposite sign of \(\beta\).

Let \(L=\partial_{tt}+\partial_{zz}-t\partial_t-z\partial_z\)
be the Gaussian generator and \(L_z=\partial_{zz}-z\partial_z\).
By \(LT_r=-r\partial_r T_r\), the value of \(LF\)
on the right graph is \(L_zG(b)+E(t)\), with
\(E(t)=-r(g-g_\infty)\). The fixed-halfspace derivative at
\(b\) is \(g_\infty\) as above.
For \(U=F-G\), put
\[
 C(u)=\beta-G(b+u),\qquad
 D(u)=L_zG(b)-L_zG(b+u).
\]
On the graph \(U=C(u)\) and \(LU=D(u)+E(t)\). Differentiation in
\(t\) along \(z=b+u(t)\) gives
\[
 (C'(u)-U_z)(u''-t u')
   =D(u)-U_{zz}+zU_z+2U_{tz}u'+(U_{zz}-C''(u))(u')^2+E(t).
                                                               \tag{CR1}
\]
All score derivatives in this formula are evaluated at the graph.

## The signed quadratic term

Choose a large longitudinal cutoff \(S\). The difference between
the right tail and its pressure halfspace contributes
\[
 W(t,z)=\int_{s>S}
  [\psi(rz-b-u(s))-\psi(rz-b)]p_\sigma(s-rt)\,ds          \tag{CR2}
\]
to \(U\). Source differences on \(s\le S\) contribute a term
\(R(t,z)\). For bounded \(z\), \(R\) and its first two
vertical derivatives are
\(O(\exp[-(rt-S)^2/(2\sigma^2)])\) by the conditional Gaussian
tail in \(s\) and the bounded vertical Gaussian derivative
moments. Here \(p_\sigma\) is the centered Gaussian density of
standard deviation \(\sigma\).

If \(m\le1\), the opposite arm and the downward axis have longitudinal
coordinates bounded above. For \(m<1\) the longitudinal coordinate
on the opposite arm tends to \(-\infty\), and for \(m=1\) it is
bounded by its asymptote. Thus the forward half-plane past
a sufficiently large \(S\) contains only the binary right graph,
and \(U=W+R\). If \(m>1\), the opposite arm instead forms a second
graph over \(s>S\), with height
\[
       w(s)=k s+O(1),\qquad k=\frac{2m}{m^2-1}>0 .
\]
Along the opposite arm, at horizontal coordinate \(-a\)
in the original frame, the new \(s\) coordinate has derivative
\((m h'(a)-1)/\sqrt{1+m^2}\), eventually positive. The right
and opposite graph heights are ordered, and no part of the
downward boundary of labels 1 and 2 enters this half-plane.
The labels in increasing vertical order are 1, 3, and 2.
Consequently \(U=W-V+R\) with
\[
 V(t,z)=T_r1_{\{s>S,\,v>w(s)\}}(t,z).                     \tag{CR3}
\]
The subtraction accounts for label 2 having pair contrast zero.

For a scalar parameter \(z\) near \(b\), write \(\xi=rz-b\) and set
\[
 k_0(z)=-\frac{r^2\psi'''(\xi)-z r\psi''(\xi)}{\psi'(\xi)},
 \qquad \mathcal B_z=-\partial_{zz}+z\partial_z-k_0(z).
\]
In the integrand of \(\mathcal B_z W\) the linear term in \(u(s)\)
cancels. Direct differentiation and Taylor expansion give
\[
 \mathcal B_z W
     =\int_{s>S}\bigl[\zeta(z)u(s)^2+O(|u(s)|^3)\bigr]
                                   p_\sigma(s-rt)\,ds ,
 \qquad
 \zeta(b)=-\frac{b r}{2(1+r)^2}\psi'((r-1)b).            \tag{CR4}
\]
The quadratic coefficient before simplifying is
\[
 \tfrac12[-r^2\psi^{(4)}+z r\psi'''-k_0(z)\psi''](\xi).
\]
The Gaussian derivative identities give its value at \(b\) as
\(-\psi'(\xi)[r^2\xi/\sigma^4+br/(2\sigma^2)]\) at
\(\xi=(r-1)b\), equal to (CR4). The Taylor remainder is uniform
for \(z\) near \(b\). Since \(b\ne0\), take \(S\) sufficiently large
that \(\mathcal B_z W\) has the sign of \(-b\) and its magnitude
bounds a positive multiple of the integral of \(u(s)^2\),
for large \(t\) and \(z=b+u(t)\).

Along the graph \(D(u)+k_0(z) C(u)\) vanishes at zero and its
derivative there is \(r\psi'((r-1)b)\). Meanwhile
\(C'(u)-U_z\to-r\psi'((r-1)b)<0\), since \(u\to0\) and the
conditional derivatives of \(W,R\), and \(V\) when present,
tend to zero. Moving the terms proportional to \(u\) or \(u'\)
in (CR1) to the left gives
\[
 u''+p(t)u'+q(t)u=f(t),\qquad
 p(t)=-t+O(1),\quad p'(t)=O(1),\quad q(t)\to1,\qquad
 f(t)=\frac{\mathcal B_z U(t,z(t))+E(t)}{C'(u(t))-U_z(t,z(t))}.
                                                               \tag{CR5}
\]
For \(p\), the added coefficient is
\(-[2U_{tz}+(U_{zz}-C'')u']/(C'-U_z)\).
The score and graph derivative bounds give \(p'=O(1)\).
For \(q\), divide the analytic function \(D+k_0C\) by \(u\),
extending continuously across zero, and use the derivative
just computed.

If \(f\) is eventually of one strict sign, then \(u\) is
eventually of one weak sign. Indeed \(v=e^{\frac12\int^t p}u\)
tends to zero and satisfies
\[
 v''+(q-p'/2-p^2/4)v=e^{\frac12\int^t p} f .
\]
The coefficient of \(v\) is negative on a tail. With positive
forcing there is no bounded positive component of \(v\) wholly
in that tail, since its positive maximum would contradict the
equation. Negative forcing is handled by replacing \(v\) by
\(-v\). Thus infinitely many changes of sign are impossible.

## When the opposite arm stays outside the forward half-plane

If \(m\le1\) and \(u\) is not eventually zero, choose a fixed
interval \([a,a']\subset(S,\infty)\) where \(u(s)^2\) is bounded
below by a positive constant, with \(a\) far enough that
\(ra/\sigma^2>C'\) in (CR0). For \(rt>a'\) the integral in
(CR4) is bounded below in magnitude by a constant times
\(p_\sigma(a-rt)\). This dominates both \(\mathcal B_z R\)
and the bound on \(E\) in (CR0). The signs in (CR4)–(CR5)
make \(f\) eventually of one strict sign. Therefore \(u\)
is eventually of one weak sign. The right tail is a bounded
isolated end approaching its pressure line from one side, so
it is exactly straight by the isolated-end theorem.
[Exact straightening of a one-sided end](THREE_CELL_ISOLATED_END_RIGIDITY.md#a-one-sided-approach-forces-exact-straightness)

If the tail was already straight the same conclusion holds.
Analyticity extends the affine formula along the whole right
arm. Reflection then gives a conical maximizing partition
about its junction. The translated-cone rigidity theorem makes
it the centered regular simplex, contradicting \(\beta\ne0\).
[Translated-cone rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

## When the opposite arm enters the forward half-plane

In the remaining case \(m>1\) of (CR0), we assume \(\beta>0\),
so \(b<0\)
and \(\mathcal B_z W\ge0\) on the tail by (CR4).
Put
\[
                 \alpha=\frac{r^2 k^2}{2\sigma^2(1+k^2)}>0 .
\]
For bounded \(z\), (CR3) is at most
\(C_\epsilon e^{-(\alpha-\epsilon)t^2}\) for every
\(\epsilon>0\). Indeed \(w(s)\ge ks-C\), so a source
point in (CR3) must satisfy \(v-ks\ge-C\). This linear
functional of the conditional Gaussian has mean
\(rz-krt\) and variance \(\sigma^2(1+k^2)\).

Applying \(\mathcal B_z\) inside \(V\) gives a polynomial
multiplier in the source normal coordinate \(v\). Uniformly
for \(z\) near \(b\), its leading coefficient is
\(-r^2/\sigma^4\), and its lower coefficients are bounded.
For all sufficiently large positive \(v\), it is negative
and bounded above by \(-c_1v^2\). On the source region
in (CR3) where \(v\) is at most this fixed threshold, both
\(s\) and \(v\) are bounded because \(s>S\) and
\(w(s)=ks+O(1)\). Its contribution is bounded in absolute
value by \(\exp[-r^2 t^2/(2\sigma^2)+O(t)]\).

For a lower bound, choose a fixed-size source rectangle
just above the line \(v=ks+C\) near the projection of
\((rt,0)\) onto \(v=ks\). Its coordinates are
\[
        s=\frac{rt}{1+k^2}+O(1),\qquad
        v=\frac{krt}{1+k^2}+O(1).
\]
For a sufficiently large fixed \(C\) the rectangle is contained
in (CR3). Its conditional density is
\(\exp[-\alpha t^2+O(t)]\), uniformly for bounded \(z\).
On it the multiplier for \(\mathcal B_z(-V)\) is positive
and at least \(c_2t^2\). Since \(\alpha<r^2/(2\sigma^2)\),
the bounded-source negative contribution, \(\mathcal B_z R\),
and \(E\) in (CR0) are negligible. Together with \(\mathcal B_z W\ge0\)
and the negative denominator in (CR5), this gives
\[
                  f(t)\le-c_\epsilon e^{-(\alpha+\epsilon)t^2}
                                 \qquad(t\ \text{large})       \tag{CR6}
\]
for each small fixed \(\epsilon>0\).
In particular \(u\) is eventually of one weak sign by
the differential equation.

Equation (CR6) gives a lower bound on the size of \(u\)
on every far interval. Integrate (CR5) against a fixed
nonnegative smooth bump of positive integral translated
to have support in \((a-1,a+1)\). By integration by parts,
the integral is a pairing with \(u\) against the test function
\(\chi''-(p\chi)'+q\chi\), bounded in magnitude by
\(C(1+a)\). It follows that
\[
           \int_{a-1}^{a+1}|u(s)|\,ds
                      \ge c'_\epsilon e^{-(\alpha+2\epsilon)a^2}
                                      \qquad(a\ \text{large}). \tag{CR7}
\]
The change from \(\epsilon\) to \(2\epsilon\) absorbs the
polynomial factor and the endpoint displacement.

Increase \(S\) to lie in the sign-definite tail. The discarded
bounded strip in \(s\) only changes \(R\) by a left Gaussian
error. For \(u\) on this tail the mean-value formula in the
stationary identity \(G(b+u)+W-V+R=G(b)\) gives the following
comparison for any small fixed \(\eta>0\) and large \(t\):
\[
 r(1+\eta)|u(t)|
 \ \ge\ (1-\eta)\int_{s>S}|u(s)|p_\sigma(s-rt)\,ds
                  -C_\eta[V(t,z(t))+|R(t,z(t))|].           \tag{CR8}
\]
The derivatives of \(\psi\) in both differences are
within \(1\pm\eta\) of \(\psi'((r-1)b)\), and the integrand
in \(W\) has one sign because \(u(s)\) does.
Taking \(a=rt\) in (CR7) bounds the integral in (CR8)
below by \(c''_\epsilon e^{-r^2(\alpha+2\epsilon)t^2}\).
Choose \(\epsilon\) so small that \(r^2(\alpha+2\epsilon)
<\alpha-\epsilon\). The error in (CR8) is then at most
\(\eta\) times the integral for large \(t\).

Choose \(\eta\) with \(r(1+\eta)<1-2\eta\), and retain
only the source interval \([\ell t,u_0 t]\) for fixed
\(0<\ell<r<u_0<1\). Inequality (CR8) bounds \(|u(t)|\)
below by more than the retained Gaussian average.
The full absorbed inequality first gives strict positivity
of \(|u(t)|\) throughout the tail, by (CR7). Choose \(R_0\)
so large that \([\ell R_0,R_0]\) is in this positive tail.
Iterate the restricted transitions from \(t>R_0\), stopping
at the first step below \(R_0\). The probability mass retained at
each step is at least \(1-Ce^{-c t^2}\) at that step.
Since the retained transitions decrease \(t\) by a fixed
factor at most \(u_0\), the product of these masses is
uniformly bounded below over all stopped paths. The
positive minimum on the stopping interval gives a uniform
positive lower bound on \(|u|\) on the tail, contrary to
\(u(t)\to0\). This excludes (CR0) in the stated slope and
pressure cases.

## A convex tail forces simplex rigidity

Suppose instead that \(h\) is convex on some tail, with no
hypothesis on \(g\). Its bounded positive derivative
has a finite positive limit \(m\). The discretely homogeneous
profile is then linear, and (TP1) gives \(h(x)=mx+c+o(1)\).
The function \(h(x)-mx\) is nonincreasing on the convex tail,
so it approaches \(c\) from above. In the rotated coordinates,
the same sign persists for \(u(t)\).

For \(m\le1\), the binary end is isolated in a forward
half-plane as in (CR2) and is exactly straight by the
one-sided end theorem. For \(m>1\), the decomposition
\(F=G+W-V+R\) is unchanged, with \(V>0\) and the Gaussian
left error \(R\). The rectangle estimate above gives
\(V\ge c_\epsilon e^{-(\alpha+\epsilon)t^2}\), so
\(V>|R|\) for all large \(t\) because
\(\alpha<r^2/(2\sigma^2)\). In the stationary identity
\(G(b+u)-G(b)=-W+V-R\) the terms from \(u(s)\ge0\) in
\(-W\) are nonnegative. The same mean-value comparison
as in (CR8), now keeping the favorable sign of \(V-R\),
gives
\[
 r(1+\eta) u(t)
             \ge (1-\eta)\int_{s>S}u(s)p_\sigma(s-rt)\,ds .
                                                               \tag{CR10}
\]
In fact the strict contribution of \(V-R\) makes \(u(t)>0\)
on the whole tail. Choose \(\eta\) small and use the same
stopped averaging argument to force a positive lower
bound on \(u\), contradicting its convergence to zero.
Thus only the straight-tail case occurs. Analytic
continuation and translated-cone rigidity give the simplex.

If \(h\) is concave on a tail, the same limit-slope argument
gives an asymptote approached from below. At \(m\le1\)
the isolated-end argument again gives the simplex. Thus
an eventually concave nonsimplex reflected maximizer would
need \(m>1\).
