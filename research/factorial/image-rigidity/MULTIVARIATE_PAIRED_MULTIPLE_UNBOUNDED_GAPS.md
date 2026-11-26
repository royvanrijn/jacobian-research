# All paired multiple Gaussian contacts with unbounded retained gaps

An asymmetric ordinary addition supplies a pure cubic operator at
the scale of an integer row split. Canceling the split's cubic at
the upper order covers the remaining Hermite rank configurations
where the lower root is double.

**Theorem.** In the
[fixed-terminal reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data),
fix a physical interior contact at \(T_*>0\), \(|\lambda_*|<T_*\)
and adjacent orders \(n,N=n+1\), \(n\ge2\). Suppose both contacting
nonreal roots have multiplicity at least two and the countable
ordinary retained core has summable curvature and unbounded
absolute gaps. Then the contact is the limit of simple transverse
interior Gaussian contacts of finite ordinary cores. The orders and
total signed gap are preserved, as are the limiting strength, bias
and relative physical root scales. Any strict total curvature cap
persists. No condition on the angular velocities is required.

## The outstanding Hermite configurations

Use the positive Hermite ranks \(h_m^-,h_m^+\) of the
[integer row splitting theorem](MULTIVARIATE_INTEGER_ROW_SPLIT_UNFOLDING.md).
That theorem gives the conclusion whenever \(h_{m_n}^+/\sqrt n>
h_{m_N}^-/\sqrt N\), in particular when \(m_n\ge3\) or \(m_N\) is
even. For \(m_n=2\) and odd \(m_N\ge11\), the rank comparison also
holds. Indeed \(\operatorname{He}_{11}'(0)=-10395\) and
\(\operatorname{He}_{11}(1)=936\), so \(h_{11}^-<1\).
Same-parity Laguerre interlacing shows \(h_{m_N}^-\le h_{11}^-\).
The case \((m_n,m_N)=(2,3)\) is covered by the
[reverse heat shift](MULTIVARIATE_REVERSE_CLOCK_ROW_REPLACEMENT.md).
It remains to prove the claim for \(m_n=2\) and
\(m=m_N\in\{5,7,9\}\) when the rank comparison fails.

For these odd \(m\), equality in the comparison cannot occur.
If it did, the squared Hermite root \((h_m^-)^2=N/n\) would
be a root of the monic integer polynomial in \(X^2\) given by
\(\operatorname{He}_m(X)/X\). This rational root is nonintegral
since \(N=n+1\) and \(n\ge2\). We therefore need only handle
\[
                        \frac1{\sqrt n}<\frac{h_m^-}{\sqrt N}.
                                                               \tag{1}
\]
The argument below for these three pairs requires only infinitely
many retained rows of absolute gap at least two.

## An ordinary pure cubic perturbation

As in the splitting theorem, take a sequence of weak rows with
positive gaps \(d\ge2\), midpoints \(M\), and \(d^2/M\to0\).
A global reversal of the columns is available to make the selected
gaps positive. Set \(\eta=d/M\) and split each row into two
almost homothetic integer rows. The local operator in \(X\)
at effective strength \(T=T_*+\eta^2y\) is
\[
                    e^{y\partial_X^2/(2j)}
                              (1+\partial_X/j)^{1/2}
                                                               \tag{2}
\]
at order \(j\), up to radial translations and an overall scalar.

Fix a small \(w>0\), take a half-integer midpoint
\(H\sim wN/\eta^3\), and add \(\lfloor w H\rfloor\) increasing
unit-gap rows with that midpoint. Put \(u=1/H\) and
\(w_H=\lfloor wH\rfloor/H\). The actual Gaussian strength is
\(T-w_H\). By the
[compensated chord formula](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients),
the log multiplier per added row begins with
\(u^2(j\psi(k/j)+\phi(k/j)/2)\) in the notation there.
With \(v=k-j/2\), its cubic coefficient is
\(-u^2 v^3/(6j^2)\). The remaining terms of this layer have degree
at most two in \(v\). At every Taylor layer \(u^q\) with \(q\ge3\),
the degree in \(v\) is at most \(q+1\), because this layer is a
finite sum of the \(q\)-th powers from the rising product.
The endpoint subtraction and compensation respect these bounds.

Now \(\lfloor wH\rfloor/(H^2\eta^3)\to1/N\).
After dilation in \(X\), lower-degree terms of the quadratic
layer vanish and the \(q\)-th layers with \(q\ge3\) have scaled
coefficients \(O_{w,j,q}(\eta^{2q-4})\), which also vanish.
As in the splitting limit, truncation of the analytic log
multipliers to a fixed sufficiently high order makes the remaining
errors on the finitely many original frequencies arbitrarily high
order in \(\eta\), even after dividing by \(\eta^{m_j}\).
The exact combined operator in the limit is consequently
\[
   e^{y\partial_X^2/(2j)-\partial_X^3/(6j^2N)}
                               (1+\partial_X/j)^{1/2}.          \tag{3}
\]
The cloud's total variance tends to \(w\), while all its relative
gaps tend to zero.

Write \(U=jX\), and make the real centering translation that
removes the linear term of \(\tfrac12\log(1+\partial_U)\).
The quadratic term of the centered operator is
\(a_j\partial_U^2\), \(a_j=jy/2-1/4\). Its cubic term at
\(j=N\) vanishes, since the cloud contributes
\(-\partial_U^3/6\) there. At \(j=n\) the root is double,
so terms of third or higher differential order annihilate the
model. Its centered polynomial is \(U^2+2a_n\), with a positive
imaginary root for \(y>1/(2n)\).

## A pair of opposite rank signs

The upper centered polynomial at \(a=a_N\) is
\[
 P_m(U,a)
  =m![z^m]\,e^{Uz+a z^2}(1+z)^{1/2}
                        e^{-z/2+z^2/4-z^3/6}.
                                                               \tag{4}
\]
Choose \(y_->1/(2n)\) close enough to \(1/(2n)\) so that
\(0<a<1/16\). At the threshold \(y=1/(2n)\),
\(a=1/(4n)\le1/36\), since \(N\ge2m\ge10\).
The exact requisite signs of (4) are as follows.

| \(m\) | \(P_m(0,a)\) | \(U_0\) | \(P_m(U_0,a)\) |
| --- | --- | --- | --- |
| 5 | \(12\) | \(1\) | \(60a^2+20a-2\) |
| 7 | \(504a+360\) | \(3/2\) | \(1260a^3+\tfrac{2835}{2}a^2-\tfrac{1953}{16}a-\tfrac{5157}{128}\) |
| 9 | \(18144a^2+25920a+15624\) | \(3\) | \(45360a^4+272160a^3+317520a^2+51840a-13050\) |

These polynomials are obtained by finite multiplication of
the series in (4). For \(0\le a\le1/16\), the middle-column
values of \(P_m(0,a)\) are positive. The upper bounds on the
last-column values obtained by discarding any negative term
of positive degree and evaluating the rest at \(a=1/16\)
are respectively \(-33/64,-35271/1024,-34826445/4096\).
Since \(P_m\) is monic of odd degree, it has at least three
real sign-changing roots. In the full imaginary-part order
at the upper cluster, rank \(1+\lceil m/2\rceil\) is therefore
zero at \(y_-\). The positive displacements in \(X\) are obtained
by dividing those in \(U\) by \(j\). The upper rank is strictly
below the positive lower rank.

At large \(y_+\), the \(\sqrt y\) scale reduces both models in
(3) to their positive-heat Hermite polynomials.
The difference of the same lower and upper ranks has leading
coefficient \(1/\sqrt n-h_m^-/\sqrt N<0\) by (1).
Use these opposite strict signs and
[ordinary regularization](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
as in the integer split proof. This gives a simple transverse
contact of a finite ordinary core after taking \(\eta\) sufficiently
small for the fixed \(w\).

The split preserves the signed degree difference, and the total
signed difference of the reference is held fixed by the residual
Gaussian data. The cloud's row center contribution tends to
\(-w/(2nN)\) by the same compensated chord formula.
At fixed \(w\), the limiting actual strength and physical bias
are thus \(T_*-w\) and \(\lambda_*-w\). Choosing
\(w<(T_*-|\lambda_*|)/2\) keeps a strict interior margin.
A diagonal choice with \(w\to0\) recovers the original strength
and bias. The transferred variance cancels the limiting added
curvature, preserving each strict total curvature bound.
Together with the existing rank and reverse-heat cases this proves
the theorem.
