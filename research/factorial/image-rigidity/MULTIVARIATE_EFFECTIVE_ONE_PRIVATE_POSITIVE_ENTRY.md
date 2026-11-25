# Effective entry at the positive one-private boundary

Use ordinary profiles with one private row \((0,1)\), other unequal
rows shared, and
\[
 \sigma=\sum_i d_i^2/a_i,\qquad
 \xi_j=\log(j^\Delta t)+\sum_i d_i\log a_i
\]
as in the [one-private comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md).

**Theorem.** Given rational \(b>0\), there are computable
positive rational \(\delta,\theta_0\) and an integer \(J\)
such that, for every profile with \(0\le\sigma\le1+\delta\),
an order-\(j\) root with \(j\ge J\) and
\(|\Im\xi_j|\le\theta_0\) must satisfy \(|\Re\xi_j|<b\).
This holds for the actual single-order moment polynomial,
uniformly in the shared rows and total degrees.

With \(b\) and \(\theta_0\) smaller if necessary, the
[effective near-endpoint theorem](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md)
therefore excludes any consecutive common root on this entire
positive-boundary sector for sufficiently large orders. There is
no lower bound on \(\sigma\).

## Coefficient ratios and the negative real fields

For \(\sigma\le5/4\), each \(r_i=d_i/a_i\) obeys
\(-5/8\le r_i\le5/4\), as in the
[bounded cubic calculation](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md).
Thus the shared quantity
\(V=\sum_i d_i^2/\min(a_i,b_i)\) is at most
\((8/3)\sigma\le10/3\).
Set
\[
 B(q)=\sum_i d_i\log(1+r_i q),\quad
 X(q)=\log(1-q)+B(q),\quad
 R(q)=(1-q)\sum_i\frac{d_i^2}{a_i+d_i q}.
\]
On the real interval \(0\le q<1\), one has
\(0\le B(q)\le10q/3\), \(X(0)=0\), and
\[
   X'(q)=\frac{R(q)-1}{1-q},\qquad
   R'(q)=-\sum_i\frac{d_i^2 b_i}{(a_i+d_iq)^2}
                     \le-\frac{2}{27}\sigma .
\]
Here \(b_i/a_i\ge3/8\) and \(1+r_iq\le9/4\).
If \(\sigma\le1\), then \(X(q)\le0\). Otherwise, put
\(\gamma=2/27\) and choose \(\delta\le\gamma/2\).
The potential increase of \(X\) is confined to
\(q\le\delta/\gamma\), where \(X'\le2\delta\).
Consequently \(X\le 2\delta^2/\gamma\).
Choose \(\delta\) small enough that this is at most \(b/4\).

Denote the exact normalized summand by \(f_{j,k}\), with
\(f_{j,0}=1\). For \(t\ne0\) the consecutive ratio is
\[
    \frac{f_{j,k+1}}{f_{j,k}}
       =\exp\bigl(\xi_j+X(k/j)+E_{j,k}\bigr),
           \qquad 0\le E_{j,k}\le V/j,\quad 0\le k<j.
                                                               \tag{2}
\]
Indeed the private factorial cancels \(k!\), so the falling
factorial contributes \(j-k\). For \(x_i=ja_i+kd_i\),
an increasing shared row has ratio
\(\prod_{s=1}^{d_i}(x_i+s)\), whose logarithm differs from
\(d_i\log x_i\) by a term between zero and \(d_i^2/(ja_i)\).
A decreasing row uses the inverse product of \(x_i-s\),
\(0\le s<|d_i|\). Its error has the same sign and is
bounded by \(d_i^2/(j\min(a_i,b_i))\) since
\(x_i-|d_i|\ge j\min(a_i,b_i)\) for \(k<j\).
The leading factors give exactly \(e^{\xi_j+X(k/j)}\).

If \(x=\Re\xi_j\le-b\), then for \(V/j\le b/4\) all such
ratios have modulus at most \(r_b=e^{-b/2}<1\).
On \(0\le q\le1/2\), \(|X(q)|\le6q\).
For \(k\le j/2\), comparing the \(k\)-fold product in (2)
to \(e^{k\xi_j}\) by telescoping gives
\[
      |f_{j,k}-e^{k\xi_j}|
               \le C(b) r_b^{k-1}(k^2+k)/j.
\]
In fact each ratio difference is at most
\(C(b)(k+1)/j\), using \(|e^u-1|\le |u|e^{\max(0,u)}\)
for real \(u=X+E_{j,k}\), and both products use the same
modulus bound \(r_b\).
The rest of the sum is bounded by its geometric tail.
Thus the normalized moment is
\((1-e^{\xi_j})^{-1}+O_b(1/j)\), with computable error,
uniformly in the imaginary field. The leading term has modulus
at least \(1/2\), so no large-order root occurs in this region.

## A single effective saddle at the positive real fields

For completeness, there is an absolute real-field cap for the
roots on the band being used. The
[bounded-private annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
gives \(|\log|G_j t||\le\log100+V/2\), with \(G_j\) the
endpoint scale. Also
\[
   \left|\log G_j-\Delta\log j-\sum_i d_i\log a_i\right|
                     \le 1+V/2+V/j .
\]
To see this, for a shared row express its endpoint shift as
\(\int_{a_i}^{b_i}\psi(1+jx)\,dx\) and bound the difference
from \(d_i\log(ja_i)\) by
\(d_i^2/(2\min(a_i,b_i))+|d_i|/(j\min(a_i,b_i))\).
Here the digamma bound
\(0\le\psi(1+y)-\log y\le1/y\) follows from log convexity
of Gamma and its unit-step difference.
The private row contributes \(\log j!/j-\log j\),
which is between \(-1\) and zero. Together these show that
one may restrict to \(|\Re\xi_j|\le13\). If \(b>13\) this
already gives the stated condition, so assume \(b\le13\).

For real \(x\in[b,13]\), let
\(\Phi(q,x)=xq+\int_0^q X(u)\,du\).
Its derivative starts positive and ends negative, and vanishes
at exactly one \(q_s\in(0,1)\).
Indeed \(R(q)\) decreases, so \(X\) first increases or
decreases, and then decreases strictly to negative infinity.
At its root \(\Phi_q=0\), the bounds on \(B(q)\) give
\[
            1-e^{-b}\le q_s\le1-e^{-13-10/3}.
\]
Choose a rational \(\eta>0\) with \(q_s,1-q_s\ge2\eta\)
uniformly by these formulas.
If \(\sigma\le1/2\) then \(R(q_s)\le1/2\).
Otherwise
\(R(q_s)\le1+\delta-q_s/27\le1-\eta/27\)
by further requiring \(\delta\le\eta/27\).
It follows that \(X'(q_s)\le-\kappa\) with
\(\kappa=\min(1/2,\eta/27)>0\).

The functions \(\Phi''=X'\) and their fixed derivatives
have explicit uniform bounds on a complex neighborhood of
\([\eta,1-\eta]\). In fact its row sum is bounded by
\((8/3)\sigma\) and the entropy derivative by the inverse
distance to one on this strip. Cauchy's estimate on a
smaller neighborhood supplies the derivative bounds.
Choose a rational \(u>0\) so small that on the complex
disk \(|q-q_s|\le4u\),
\(|X'(q)-X'(q_s)|\le\kappa/4\).
At the real points \(q_s\pm u/2\), the real phase is
below its maximum by at least \(\gamma_2=\kappa u^2/32\),
by integrating the negative Hessian.
Uniqueness and the sign of the real derivative extend
this gap to the whole real segment outside that interval.

For \(\xi=x+i\theta\) with \(|\theta|\) sufficiently small,
Rouché on \(|q-q_s|=u/4\) gives the unique nearby complex
saddle \(q_\theta\), and
\(|q_\theta-q_s|\le2|\theta|/\kappa\).
Indeed \(\Phi_q(q,\xi)=\Phi_q(q,x)+i\theta\) and the
uniform derivative bound compares \(\Phi_q(q,x)\) with
\(X'(q_s)(q-q_s)\).
This also gives \(\Re(-\Phi_{qq}(q_\theta,\xi))\ge\kappa/2\).
The real saddle phase changes by at most \(C\theta^2\)
from \(\Phi(q_s,x)\), since \(\Re\Phi(q_s,x+i\theta)
=\Phi(q_s,x)\) and the change to \(q_\theta\) has size
\(O(|\theta||q_\theta-q_s|+|q_\theta-q_s|^2)\).

## Uniform actual lattice estimate

Keep the real lattice interval between half-grid endpoints
within \(1/j\) of \(q_s-u,q_s+u\).
The omitted sum divided by \(j\) is exponentially small
relative to \(e^{j\Phi(q_s,x)}\) by the real phase gap.
The summand bound \(Cj e^{j\Re\Phi(q,\xi)}\) on the whole
integer grid follows from the falling-factorial Riemann estimate
and shared Stirling bounds in the
[effective cubic expansion](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md#uniform-real-lattice-expansion).
Those estimates apply on this band without any small-field
assumption for this upper bound, since the field enters the
summand as the exact factor \(e^{jq\xi}\).

On the retained interval the analytic summand is
\(A_j(q)e^{j\Phi(q,\xi)}\), with the same uniform Stirling
bounds as in the
[effective cubic contour calculation](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#effective-lattice-representation-on-a-fixed-small-rectangle).
Choose a fixed rational contour height \(\nu>0\) so small
that the phase increases by less than \(\gamma_2/4\)
on the endpoint walls at \(\theta=0\), and by less than
\(\pi\nu/4\) on the horizontal sides.
These choices are computable from the uniform second-derivative
bound and the real values of the first derivative.
Now require \(|\theta|\nu\le\min(\gamma_2/4,\pi\nu/4)\).
The cotangent kernel at half-grid endpoint walls is bounded;
on horizontal sides its nonconstant part decays like
\(e^{-2\pi j\nu}\). The constant parts give the real
integral as in the cubic contour calculation.
Thus the normalized lattice sum divided by \(j\) equals
the retained integral up to \(Cj^C\exp(j\Phi(q_s,x)-\gamma_3 j)\)
with a computable \(\gamma_3>0\).

Decrease \(\theta_0\) further so that \(q_\theta\) stays well
inside the retained interval's complex neighborhood,
\(|\Im q_\theta|\le\nu\), and its real phase loss
\(C\theta^2\) is less than \(\gamma_3/4\).
Shift the retained integral horizontally to the line through
\(q_\theta\), connecting at its two ends. Those short endpoint
connectors are exponentially small by the same real gap
and derivative bounds. On the horizontal segment the
real part of the second derivative is \(\le-\kappa/2\),
so its maximum is the continued saddle and has a uniform
quadratic loss. Taylor expansion near the saddle gives
\[
     \frac1j\sum_{k=0}^j f_{j,k}
      = A_0(q_\theta)e^{j\Phi(q_\theta,\xi)}
             \sqrt{\frac{2\pi}{-j\Phi_{qq}(q_\theta,\xi)}}
                 (1+O(j^{-1/2})).
\]
All constants are computable. In the local Gaussian estimate
one bounds the cubic remainder by a uniform constant times
the cube of distance, integrates the resulting error against
a fixed fraction of the quadratic majorant, and uses the
strict loss on the complement. The explicit Stirling
amplitude \(A_0(q)=
(1-q)^{-1/2}\exp(\frac12\sum_i\log(1+r_iq))\)
is uniformly bounded away from zero on the saddle disks.
The displayed nonzero leading term dominates the error
above a computable \(J\).
This completes the positive real field exclusion.

Taking \(b<\xi_*/2\), \(\theta_0<\xi_*/2\), and
\(\delta\le\delta_*\) from the effective endpoint theorem
forces both fields of a common root into its regime.
Their imaginary parts are equal at the two orders, because
the exact physical drift is real. This gives the claimed
effective consecutive exclusion near the positive real axis.
