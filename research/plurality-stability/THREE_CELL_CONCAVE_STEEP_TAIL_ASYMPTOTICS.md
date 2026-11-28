# Gaussian rates on an eventually concave steep reflected tail

Fix \(0<r<1\). Let a balanced planar three-cell Gaussian maximizer
have the central reflection interchanging labels 1 and 2. Use
the [reflected graph coordinates](THREE_CELL_REFLECTION_REDUCTION.md)
with right height \(h(x)\), and suppose \(h\) is eventually concave
with \(\lim_{x\to\infty}h'(x)=m>1\). Set
\[
   \sigma^2=1-r^2,\qquad k=\frac{2m}{m^2-1},\qquad
   \lambda=\frac{k^2}{2\sigma^2},\qquad
   \alpha=\frac{r^2 k^2}{2\sigma^2(1+k^2)} .
\]
In orthonormal coordinates \((t,z)\) running along and normal to
the right affine asymptote, the graph is \(z=b-q(t)\) on a tail,
where \(q>0\) is nonincreasing. Then
\[
   \lim_{t\to\infty}t^{-2}\log q(t)=-\lambda .
                                                               \tag{CS1}
\]
Let \(g(t)\) be the correlation derivative of the pair score
\(T_\rho(1_{A_3}-1_{A_1})\) at \(\rho=r\), with the partition
and the observation \((t,b-q(t))\) held fixed. Let \(g_\infty\)
be the corresponding derivative for the binary pressure halfspace
at \((t,b)\). The score difference \(g(t)-g_\infty\) is positive
for all large \(t\) and has Gaussian rate
\[
             \lim_{t\to\infty}t^{-2}\log(g(t)-g_\infty)
                       =-\alpha .                          \tag{CS2}
\]
In fact, the leading coefficient relative to the conditional
opposite-arm contribution is given in (CS5) below. These are
necessary conditions on the eventually concave steep case.
That case also requires \(k<r\), or equivalently
\(m>(1+\sqrt{1+r^2})/r\).

## Positive contributions to the stationary score

The finite slope limit makes the phase profile \(H\) linear,
and (TP1) gives \(h(x)=mx+c+o(1)\). Eventual concavity implies
that \(h(x)-mx\) is nondecreasing to \(c\). In the rotated
coordinates the normal intercept is \(b=c/\sqrt{1+m^2}\),
and \(q\ge0\) is nonincreasing to zero because the longitudinal
coordinate increases along the graph. If \(q=0\) on a tail,
analyticity extends that affine arm to the junction, and
translated-cone rigidity would make \(m=1/\sqrt3\).
Thus \(q>0\) on the tail.
[Pressure correction, TP1](THREE_CELL_REFLECTED_TAIL_PHASE_REGULARITY.md)
[Translated-cone rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

For a large fixed \(S\), the opposite arm on \(s>S\) is a graph
with height \(w(s)=ks+O(1)\). The three labels in order of
increasing normal coordinate are 1, 3, and 2. Indeed the
longitudinal coordinate of the opposite arm at original
horizontal coordinate \(-a\) is \((-a+mh(a))/\sqrt{1+m^2}\).
It has an eventually positive derivative and tends to infinity,
and the normal coordinate is \((ma+h(a))/\sqrt{1+m^2}\).
The downward axis has bounded longitudinal coordinate above.

Write \(\psi(u)=2\Phi(u/\sigma)-1\) and \(G(z)=\psi(rz-b)\).
At observation \((t,z)\), write the conditional source density
as \(p_\sigma(s-rt)p_\sigma(v-rz)\), with \(p_\sigma\) the
centered normal density of standard deviation \(\sigma\).
The actual score has the decomposition \(F=G+W-V+R\), where
\[
 \begin{split}
 W(t,z)&=2\int_{s>S}\int_{b-q(s)}^b
                      p_\sigma(s-rt)p_\sigma(v-rz)\,dv\,ds,\\
 V(t,z)&=\int_{s>S}\int_{w(s)}^\infty
                      p_\sigma(s-rt)p_\sigma(v-rz)\,dv\,ds .
 \end{split}                                                \tag{CS3}
\]
The score remainder \(R\) comes from \(s\le S\) and is bounded
by \(C\exp[-r^2t^2/(2\sigma^2)+O(t)]\) for bounded \(z\).
Here and below \(z=b-q(t)\) unless stated otherwise.
The same bound holds for the fixed-observation correlation
derivative of \(R\), since differentiating the conditional kernel
inserts a polynomial and leaves the horizontal Gaussian tail.

The opposite term satisfies
\(\log V(t,z)=-\alpha t^2+O(t)\).
For the upper bound use \(w(s)\ge ks-C\) and the Gaussian
probability of \(v-ks\ge-C\). For the lower bound use a
fixed-size rectangle above \(v=ks+C\) at
\[
          s=\frac{rt}{1+k^2}+O(1),\qquad
          v=\frac{krt}{1+k^2}+O(1).
\]
In particular \(R=o(V)\) since \(\alpha<r^2/(2\sigma^2)\).
The stationary identity at \(z=b-q(t)\) gives
\[
                 W=V-R+G(b)-G(b-q(t)).                   \tag{CS4}
\]
Let \(d_0=\psi'((r-1)b)>0\). For any small \(\eta>0\),
increasing \(S\) makes \(W\) bounded below by
\((1-\eta)d_0\int_{s>S}q(s)p_\sigma(s-rt)\,ds\) and above
by the analogous bound with \(1+\eta\).
This follows from the vertical density in (CS3) and the
smallness of \(q(s)\) and \(q(t)\). Also
\(G(b)-G(b-q(t))\le r(1+\eta)d_0 q(t)\).
By monotonicity, \(q(s)\ge q(t)\) on \(S<s\le t\),
whose conditional horizontal probability tends to one.
Taking \(\eta\) small compared with \(1-r\) in (CS4)
therefore gives \(q(t)\le C V\) and
\(c V\le W\le C V\) for positive constants \(c,C\).

## The height rate and the source locations

For large \(s\), keep the interval \([s-1,s]\) in the
integral for \(W\). Its heights are at least \(q(s)\).
Choosing \(t=(1+k^2)s/r\), the horizontal density on this interval
is bounded below by \(c\exp[-(k^2s+1)^2/(2\sigma^2)]\).
The upper bound \(W\le CV\) gives
\[
                    q(s)\le C e^{-\lambda s^2+C's}.       \tag{CS6}
\]
The coefficient is
\(\alpha(1+k^2)^2/r^2-k^4/(2\sigma^2)=\lambda\).
Since \(\lambda>\alpha\), we have \(q(t)=o(V)\), and
(CS4) sharpens to \(W/V\to1\).

Normalize the positive integrands of \(W,V\) to probability
densities on the conditional source points. Under the \(W\) density,
\((s/t,v/t)\) concentrates at \((a_*,0)\), where
\(a_*=r/(1+k^2)\). The vertical coordinate is bounded, and the
horizontal upper rate from the last display is
\[
       J_\epsilon(a)=\frac{(a-r)^2}{2\sigma^2}
                              +(\lambda-\epsilon)a^2 .
\]
At \(\epsilon=0\) the unique minimizer over \(a\ge0\) is
\(a_*\) and its value is \(\alpha\). Outside any neighborhood
of \(a_*\), that value has a strict positive gap for
sufficiently small \(\epsilon>0\). This upper estimate, integrated
outside the neighborhood, is negligible compared with \(W\),
which has rate \(\alpha\). The integration costs only a
polynomial factor since the quadratic is coercive.

For the normalized \(V\) density, \((s/t,v/t)\) concentrates
at \((a_*,ka_*)\). The objective
\(((a-r)^2+b_1^2)/(2\sigma^2)\) on
\(\{a\ge0,\ b_1\ge ka\}\) has that unique minimizer with
value \(\alpha\). The \(O(1)\) errors in \(w(s)\) and in
the normal conditional mean vanish after rescaling.
Outside a neighborhood of the minimizer the same strict
quadratic gap and Gaussian integration give concentration.
For both densities the second moments converge as well.
On the complement of a sufficiently large bounded set of
scaled source coordinates, the Gaussian density times any
fixed quadratic factor has a rate larger than \(\alpha\),
and the denominators have rate \(\alpha\).

For the lower bound on \(q\), fix \(0<\varepsilon<\delta\) with
\(\varepsilon<k^2(1+\delta)\) and let
\(t=(1+k^2)(1+\delta)s/r\). The \(W\) saddle then has horizontal
coordinate \((1+\delta)s\). Its concentration puts
\(1-o(1)\) of the mass on the interval
\([(1+\delta-\varepsilon)s,(1+\delta+\varepsilon)s]\).
On this interval the source height deficit is at most \(q(s)\)
by monotonicity, and the horizontal Gaussian density is bounded
above by its value at the nearer endpoint to \(rt\). It follows
that
\[
 \liminf_{s\to\infty}s^{-2}\log q(s)
 \ge -\frac{\alpha(1+k^2)^2(1+\delta)^2}{r^2}
           +\frac{[k^2(1+\delta)-\varepsilon]^2}{2\sigma^2}.
\]
Let \(\varepsilon,\delta\) tend to zero. The resulting lower
rate is \(-\lambda\), proving (CS1).

## The correlation derivative

At fixed observation \(z=b-q(t)\), differentiation of the
conditional density with respect to the correlation multiplies it
by
\[
 M(s,v;t,z)=\frac{2r}{\sigma^2}+
       \frac{(s-rt)t+(v-rz)z}{\sigma^2}
       -\frac{r[(s-rt)^2+(v-rz)^2]}{\sigma^4}.
\]
For scaled sources \((s/t,v/t)\to(a,b_1)\),
\(M/t^2\) tends to
\((a-r)/\sigma^2-r[(a-r)^2+b_1^2]/\sigma^4\).
The second moment bounds justify taking the expectations under
the two normalized source densities. The limiting \(W\)
expectation minus the limiting \(V\) expectation is
\(rk^2a_*^2/\sigma^4\).

All derivative terms from \(R\) are negligible compared with
\(t^2 V\), by their faster rate. The difference between the fixed
binary halfspace derivative evaluated at \(z=b-q(t)\) and at
\(z=b\) is \(O(q(t))\), again negligible. Using \(W/V\to1\),
we obtain
\[
    \lim_{t\to\infty}\frac{g(t)-g_\infty}{t^2 V(t,b-q(t))}
          =\frac{r^3 k^2}{\sigma^4(1+k^2)^2}>0 .           \tag{CS5}
\]
This proves (CS2).

## Exclusion at and above the critical opposite slope

Suppose \(k\ge r\). Then
\(\lambda=k^2/(2\sigma^2)\ge r^2/(2\sigma^2)\).
Let \(D(s)\) be the vertical projection, at the fixed observation
height \(b\), of the full difference between the true pair contrast
and the binary contrast with dividing line \(v=b\):
\[
 D(s)=\int_{\mathbb R}
   [1_{A_3}(s,v)-1_{A_1}(s,v)-\operatorname{sgn}(v-b)]
                                    p_\sigma(v-rb)\,dv .
\]
This bounded function uses the true partition throughout the plane.
For \(P_\sigma\) denoting convolution with \(p_\sigma\), its
horizontal transform is
\[
        P_\sigma D(rt)=F(t,b)-G(b)=F(t,b)-F(t,b-q(t)).
\]
The vertical derivative of the smoothed pair contrast is uniformly
bounded by the \(L^1\) norm of the Gaussian derivative. Thus (CS6)
gives \(|P_\sigma D(x)|\le C\exp[-x^2/(2\sigma^2)+C''x]\)
for sufficiently large positive \(x\). Translate by a sufficiently
large constant \(L\), and rescale to variance one. Then
\(u\mapsto (P_\sigma D)(L+\sigma u)\) is the variance-one Gaussian
convolution of \(u\mapsto D(L+\sigma u)\), and is bounded in
magnitude by \(C e^{-u^2/2}\) for positive \(u\).
The critical one-sided Gaussian deconvolution theorem therefore
gives \(\operatorname{supp}D\subseteq(-\infty,L]\).
[Critical one-sided Gaussian deconvolution](GAUSSIAN_EXTERIOR_LINEAR_RIGIDITY.md#one-sided-gaussian-deconvolution)

On the far interval where the two arm graphs have the geometry
of (CS3), the identity \(D=0\) reads
\[
  Q(s)=P(s),\qquad
  Q(s)=\psi((r-1)b+q(s))-\psi((r-1)b),\quad
  P(s)=\bar\Phi((w(s)-rb)/\sigma).
                                                               \tag{CS7}
\]
It holds pointwise there by analyticity. Let \(t_0=mh(0)/\sqrt{1+m^2}\)
be the longitudinal coordinate of the junction. The right graph
is analytic over \(t>t_0\), so \(q(t)\), now extended as \(b\)
minus its normal height on the whole right graph, and \(Q(t)\)
are analytic throughout that interval.

On the opposite arm consider the maximal interval \(a>a_0\)
on which its longitudinal coordinate
\(s(a)=(-a+mh(a))/\sqrt{1+m^2}\) has positive derivative.
Here maximal means the component of positive derivative that
contains the tail. Its image is \((s_{\min},\infty)\), where
\(s_{\min}=s(a_0)\), and \(w(s)\) is analytic there.
By the identity theorem, (CS7) holds for all
\(s>\max(t_0,s_{\min})\). If \(s_{\min}>t_0\), then \(a_0>0\)
since \(s(0)=t_0\), and \(s'(a_0)=0\). The normal coordinate
\((ma+h(a))/\sqrt{1+m^2}\) has derivative bounded away from zero
near \(a_0\). Consequently \(w'(s)\to+\infty\) as \(s\) decreases
to \(s_{\min}\). The value \(w(s_{\min})\) is finite, so
\(P'(s)\to-\infty\), contradicting analyticity of \(Q\) at
\(s_{\min}>t_0\).

It follows that (CS7) holds on the whole interval \(s>t_0\).
There \(w'(s)>0\), so \(P'(s)<0\), and \(\psi'>0\) implies
\(q'(s)<0\). Returning to the original graph gives \(h'(x)>m\)
for every \(x>0\). In particular \(h(x)>h(0)+x/\sqrt3\)
throughout \(x>0\). This contradicts the required crossing of
that ray by a nonsimplex reflected maximizer. Thus \(k<r\).
[Crossing of the simplex ray at the junction](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#the-graph-crosses-the-simplex-ray-through-its-junction)
