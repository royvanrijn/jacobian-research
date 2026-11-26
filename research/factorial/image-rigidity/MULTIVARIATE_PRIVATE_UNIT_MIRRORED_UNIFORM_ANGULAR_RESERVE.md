# Uniform angular reserve and stability near mirrored private-unit profiles

An ordinary mirrored private-unit profile has one private row
\((1,0)\), one private row \((0,1)\), and any number of shared
mirrored pairs \((\ell,\ell+d),(\ell+d,\ell)\), with
\(\ell,d\ge1\) integers. Equal rows can be omitted.
Write \(P_m\) for the endpoint-normalized moment polynomial.
The same reference arguments below allow a common nonnegative
Gaussian coefficient strength \(\tau\), with multiplier
\(e^{-\tau k(m-k)/(2m)}\).
The estimates also apply to the mirrored pairs of real exponents
and equal positive Gamma shapes in the
[matched folded-law proof](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md).
That broader reference class is used only for a comparison
with ordinary profiles.

**Theorem.** For any such reference, at consecutive orders
\(n,N=n+1\), the circular distance between any two arguments
from the different orders is strictly greater than
\(1/(10^4 N^4)\) and \(1/(150000000N^2)\).
At every root of \(P_m\), \(m\ge2\),
one has \(\lvert P_m'\rvert\ge1/4\).
Both bounds are independent of the number and exponents of
shared rows, and of the reference curvature.

For a consequence, let \(J\ge4\) be any integer at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If an ordinary profile with one private unit row in
each orientation and only shared rows thereafter has
endpoint-normalized log coefficients within
\[
                         \frac{1}{10^7 J^5}                     \tag{1}
\]
of one such reference at every order \(2\le m<J\), then the
ordinary profile has disjoint consecutive arguments in every
degree. In particular one may adjoin arbitrary ordinary shared
rows to any mirrored private-unit profile, preserving the
all-degree conclusion whenever their summed
\(V_{\rm tail}=\sum_i (b_i-a_i)^2/\min(a_i,b_i)\) is at most
\(8/(10^7 J^6)\). No restriction is placed on the curvature
of the original mirrored core. The numerical value of \(J\)
is unevaluated.

A second log coefficient tolerance is
\(1/(2\cdot10^{11}J^3)\), which suffices independently of (1).
For a shared tail this gives \(V_{
m tail}\le4/(10^{11}J^4)\).
The two bounds are useful at different values of \(J\).

## The resolvent and its phase slopes

Use the
[mirrored private-unit comparison](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md).
For \(0<\theta<\pi\), with \(c=\cos\theta\), it gives
\[
 K_m(\theta)=A_m-B_m e^{-i\theta},\qquad
 P_m(e^{i\theta})=K_m(\theta)+e^{im\theta}\overline{K_m(\theta)} .
\]
Here \(A_m=\mathbb E[1/D]\), \(B_m=\mathbb E[T/D]\),
\(D=1-2cT+T^2\), and the expectation is under the
normalized lower endpoint tilt at degree \(m\).
The midpoint expression has \(R_m=A_m/B_m\ge1\).
With \(r_m=R_m-1\), write
\[
         \phi_m(\theta)=\arg(R_m-e^{-i\theta}),\qquad
         \psi_m(\theta)=m\theta/2-\phi_m(\theta).
\]
The positive roots lie in the angle boxes of that proof
and satisfy \(\psi_m(\theta_{m,j})=(j-1/2)\pi\).
The geometric case has \(R_m=1\), and has its positive
angles at the right endpoints of the boxes.

The expression for \(r_m\) is a quotient of two positive
integrals with the same \(\theta\)-dependent kernel
\(1/(\cosh s-\cos\theta)\). Its numerator uses
\(\cosh((m/2+1)s)-\cosh(ms/2)\) and denominator uses
\(\cosh(ms/2)\). Their ratio is increasing for \(s\ge0\).
Differentiating the common kernel, or comparing the
two reweighted expectations, shows that
\[
 0\le r_m'(\theta)\le r_m(\theta)\cot(\theta/2).
                                                               \tag{2}
\]
The logarithmic derivative of the kernel lies between
\(-\cot(\theta/2)\) and zero. Its expectation in the
numerator law is at least its expectation in the
denominator law. Subtracting yields (2). This also
holds for the geometric case \(r_m=0\).

Write \(\delta=1-\cos\theta\) and \(h=\sin(\theta/2)\).
At fixed \(R_m\), the derivative of \(\phi_m\) in \(\theta\)
has absolute value at most
\(\lvert R_m-e^{-i\theta}\rvert^{-1}\le1/(2h)\).
The term from \(R_m'\) is nonpositive and its magnitude is,
by (2), at most
\[
  \frac{r_m(1+c)}{r_m^2+2R_m\delta}
                         \le\frac1{\sqrt{2\delta}}=\frac1{2h},
\]
using \(x/(x^2+2\delta)\le1/(2\sqrt{2\delta})\).
Therefore
\[
                  |\psi_m'(\theta)|\le 3m/2
                            \quad(\theta\ge\pi/m).             \tag{3}
\]
For a lower bound, the fixed-\(R_m\) contribution is
\((R_m c-1)/[(R_m-c)^2+\sin^2\theta]\). If
\(\pi/m\le\theta\le\pi/2\), this is at most
\(1/(2\sin\theta)\le m/4\), by maximizing
\(x/(x^2+\sin^2\theta)\) in \(x\ge0\).
It is nonpositive for \(\pi/2\le\theta\le\pi\).
Hence \(\psi_m'\ge m/4\) in the first range and
\(\psi_m'\ge m/2\) in the second.

The moment \(\mathbb E T\) is the first normalized coefficient,
at most one by log convexity and the endpoint normalization.
Since \(D\le(1+T)^2\), Jensen's inequality gives
\(A_m\ge(1+\mathbb E T)^{-2}\ge1/4\).
Also \(|K_m|\ge A_m\sin\theta\) for \(0<\theta\le\pi/2\)
and \(|K_m|\ge A_m\) for \(\pi/2\le\theta\le\pi\).
At a root,
\(\lvert P_m'\rvert=2|K_m|\psi_m'\).
The angle boxes give \(\sin\theta\ge2/m\) in the first range.
The bounds above yield \(|P_m'|\ge1/4\) there and at
all other roots, including the odd structural root at
\(-1\) by continuity and the second range.

## Quantifying the adjacent phase margin

Let \(n\ge3\), \(N=n+1\), and fix a positive nonreal lower
root angle \(\theta\). For \(r=R_n(\theta)-1\), define the
comparison value
\[
          R_*=\cos\theta+\frac{(1+\cos\theta)r}{r+2}.
\]
The phase at \(R_*\), defined by the argument of
\(R_*-e^{-i\theta}\), is exactly \(\theta/2+\phi_n(\theta)\).
All these arguments are in \([0,\pi/2]\).

For \(\pi/n\le\theta\le\pi/2\), put \(t=n/N\),
\(y=\operatorname{arcosh}(1+r)\),
\(p=2+h+y\), and \(b_0=R_*-1=2(cr-\delta)/(r+2)\).
The matched phase proof gives \(R_N-1\ge r t^p\).
We sharpen its strict scalar comparison to
\[
                         r t^p-b_0\ge\frac1{1200 n^2}.         \tag{4}
\]
For \(n\ge4\) and \(b_0>0\), the scalar proof in that note
gives
\[
                         \log(r t^p/b_0)\ge\frac1{50n}.
                                                               \tag{5}
\]
Here are the margins in that inference. The logarithmic
margin is increasing in \(\theta\), so one can evaluate at
\(\pi/n\) where \(h<2/5\). For \(0<r\le5\) its lower
bound multiplied by \(n\) is
\(2\sqrt{46(r+1)/(5(r+4))}-(12/5+y)\).
The Bernstein polynomial in the cited scalar proof
has every coefficient at least \(7853/2800>14/5\).
The difference of the squares of the two displayed
terms is therefore greater than \((14/5)/(r+4)\ge14/45\).
Their sum is less than \(11\), giving more than \(1/50\).
For \(r\ge5\) the lower bound
\(4\log(1+r/2)-(12/5+y)\) starts greater than \(1/10\)
and increases. If \(b_0\) ceases to be positive,
inequality (5) is replaced by the direct positive bound
\(r t^p-b_0\ge r t^p\).

For \(0<r\le5\) and \(r\le n\delta/12\), use the identity
\[
 r t^p-b_0=\frac{r}{1+r/2}
  \left[\exp\big(\log(1+r/2)-p\log(1+1/n)\big)-1
                 +\delta(1+1/r)\right] .
\]
Here \(p<6\), and \(e^x-1\ge x\) bounds its bracket
below by \(\delta/(2r)\). Thus the margin is at least
\(\delta/7\). This also bounds the margin at \(r=0\)
by continuity. For \(n\delta/12<r\le5\), one has
\(r t^p>r/4\). If \(b_0>0\), (5) bounds the margin by
at least \(r t^p/(100n)>\delta/4800\); if \(b_0\le0\)
the bound is even larger. Since \(n^2\delta>23/5\),
both prove (4). For \(r\ge5\) the function \(r t^p\)
is increasing in \(r\), since its logarithmic derivative is
\(1/r-\log(1+1/n)/\sqrt{r(r+2)}>0\).
Its value at five exceeds one, so (5) again proves (4).

For \(n=3\) use \(\delta\ge1/2\). If \(r\le1/4\), then
\(b_0\le-\delta/(r+2)\). For \(1/4\le r\le8\),
the scalar bounds in the cited proof give
\(r t^p>r(3/4)^6\) and \(b_0<r/6\);
the margin is greater than
\(((3/4)^6-1/6)/4>1/(1200\cdot9)\).
For \(r\ge8\), the lower bound \(r t^p\) is increasing
and exceeds \(8(3/4)^6\), whereas \(b_0\le1\).
This completes (4).

If instead \(\theta\ge\pi/2\) and \(r\le5\), then
\(1-R_*=2(1-c(1+r))/(r+2)\ge2/7\).
If \(r>5\), (2) implies
\(R_n(\theta)\le2R_n(\pi/2)\), since the integral of
\(\cot(\theta/2)\) from \(\pi/2\) to \(\pi\) is \(\log2\).
Thus \(R_n(\pi/2)-1>2\).
The matched comparison there, and monotonicity of
\(r t^p\) in \(r\), give
\(R_N(\theta)\ge R_N(\pi/2)>1+2(3/4)^5>1+1/4\).
In both cases \(R_N(\theta)-R_*>1/4\).

In all cases the phase difference at the lower angle is
\[
 \psi_N(\theta)-\psi_n(\theta)
 =\int_{R_*}^{R_N(\theta)}
     \frac{\sin\theta\,dv}{(v-\cos\theta)^2+\sin^2\theta}.
\]
Restrict to the first \(1/(1200n^2)\) of the integration
interval. Since \(0\le R_*-\cos\theta\le2\), the
denominator is at most ten there. The angle boxes give
\(\sin\theta\ge2/N\). Thus the phase difference is at
least \(1/(6000 n^2 N)\).
It is also at most \(\pi/2\), since
\(\phi_n\le(\pi-\theta)/2\) and \(\phi_N\ge0\).
By strict interlacing, the lower angle is between the
upper roots with the same and next labels when both
exist. By (3) its distance to either is greater than
\(1/(10^4N^4)\). The angle boxes themselves give larger
crossing distances at zero and at the structural root
\(-1\), and reflection treats the lower half-plane.

At \(n=2\), the lower nonreal root has
\(\theta\ge\pi/2\) by its angle box. The same large-angle
comparison applies. The lower matched resolvent bound
at \(\pi/2\) uses \(t=2/3\) and
\(2(2/3)^5>1/4\) in place of \(2(3/4)^5>1/4\).
Thus \(R_N-R_*>1/4\) there and the same phase integral
gives the stated gap. At \(n=1\) the quadratic angles
are at least \(\pi/3\) from \(\pi\).
This finishes the all-degree reference gap.

## A quadratic lower gap in the order

For \(n\ge120\) and \(\pi/n\le\theta\le\pi/2\), write
\(s=\sin\theta\) and \(l=r(n/N)^p\) as above. The tangent
of the difference of the two arctangents gives
\[
 \frac{\phi_N-\phi_n}{\theta}
 \le \frac{r(1-t^p)}{r^2t^p+s^2}
 \le \min\left(\frac{p\, t^{-p/2}}{2n s},
                  \frac{e^{p\log(1+1/n)}-1}{r}\right) .
                                                               \tag{6}
\]
This uses \(R_N-1\ge l\), \(s\le\theta\), and
\(1-e^{-p\log(1+1/n)}\le p/n\).
For \(r=0\) the phase difference is nonpositive.

If \(0<r\le1/10\) and \(\theta\le\pi/6\), then
\(y<9/20\), \(h<11/42\), \(p<11/4\),
\(t^{-p/2}<51/50\), and \(n s>31/10\).
For the last two inequalities use
\(e^x\le(1-x)^{-1}\) for \(0\le x<1\) and
the cubic sine lower bound at \(\pi/n\).
The first bound in (6) is below \(47/100\).
For the same \(r\) and \(\theta\ge\pi/6\),
use \(p<7/2\), \(t^{-p/2}<51/50\) and \(2ns\ge n\),
which gives an even smaller bound.
For \(r\ge1/10\) the second bound in (6) decreases with
\(r\). In fact \(r\,\partial_r p\le1\) and, for
\(L=\log(1+1/n)\le1/120\),
\(1-e^{-Lp}\ge1-e^{-2L}\ge L\).
These give the derivative sign after multiplying by
\(r^2 e^{-Lp}\).
At \(r=1/10\), \(p<7/2\) and the bound is at most
\(70/233<47/100\), again by \(e^x\le(1-x)^{-1}\).
Thus \(\psi_N-\psi_n>3\theta/100>1/(20N)\)
on the whole first range at these orders.

For \(\pi/2\le\theta\le3\pi/4\), the bound
\(R_N-R_*>1/4\) above and the same phase integral on
\([R_*,R_*+1/4]\) give a margin greater than
\(\sin\theta/40>1/100\).
For \(\theta\ge3\pi/4\) the margin is at least
\(\theta-\pi/2\), using
\(\phi_N\le(\pi-\theta)/2\).
In particular it is greater than \(1/(20N)\)
in both ranges for \(n\ge120\).
Equation (3) now gives cross-degree gap greater than
\(1/(30N^2)\) at these orders. When \(N\le120\) use
the already proved \(1/(10^4 N^4)\), which is at least
\(1/(144000000N^2)\). This proves the second gap
in the theorem.

## Stability and a matched-slot application

All normalized coefficients of both profiles lie in \((0,1]\)
after the private factorial cancellation. Hence the log
difference bound (1) is also an absolute coefficient
difference bound. For \(m<J\), put \(\rho=1/(2\cdot10^5J^4)\)
and draw the disk of radius \(\rho\) about each reference
root. The angle boxes imply a same-degree angular mesh of
at least \(\pi/m\), so the disks are disjoint.
In \(|z|\le1+\rho\), the coefficient bound gives
\(\lvert P_m''(z)\rvert\le3m^3\), since
\((1+\rho)^m<3\). By the derivative lower bound and
Taylor's formula, on each disk boundary
\[
                  |P_m(z)|\ge\rho/4-3m^3\rho^2/2>\rho/8 .
\]
On that boundary the absolute evaluation change is at
most \(3J/(10^7J^5)<\rho/8\). Rouché assigns one
actual root to each reference disk. The angular error
is at most \(2\rho\), so the two orders cannot share
an argument in view of the reference gap.
Terminal two follows from real positivity in even
degree, and \(N\ge J\) from the absolute theorem.

For the second tolerance, use disks of radius
\(\rho=1/(4\cdot10^{9}J^2)\) at the reference roots.
From the angular mesh of at least \(\pi/m\), the sum
of the reciprocal distances from one root to the others
is at most \(m\sum_{j=1}^{\lfloor m/2\rfloor}1/j\le m^2\).
Indeed in either cyclic direction the \(j\)-th root has
chord distance at least \(2j/m\) when its corresponding
angular separation is at most \(\pi\).
The same bound follows directly with the smaller cyclic
index distance if the other direction is shorter.
For \(|w|=\rho\), factor the polynomial at a reference
root \(\zeta\) to obtain
\[
 |P_m(\zeta+w)|
 \ge \rho |P_m'(\zeta)|
   \left(1-\rho\sum_{\xi\ne\zeta}|\zeta-\xi|^{-1}\right)
 > \rho/8 .
\]
The product bound follows by multiplying factors
\(1-\rho/|\zeta-\xi|\), each in \([0,1]\).
The coefficient error with the second tolerance costs
at most \(3J/(2\cdot10^{11}J^3)<\rho/8\).
As before Rouché pairs the roots. Their two angular
errors total at most \(4\rho\), strictly below
\(1/(150000000J^2)\).

Appending a shared tail changes the absolute log
coefficients by at most \(m V_{\rm tail}/8\) by the
Gamma chord Green bound, proving the stated consequence.
In particular shared absolute curvature at most
\(4/(10^7J^6)\) suffices, since such a tail has every
relative shared gap at most one and \(V_{\rm tail}\)
at most twice its absolute curvature.

More generally, compare the shared Gauss slots of an
arbitrary original private-unit profile with those of a
mirrored reference, matching only slots of the same
orientation. A slot \((u,x)\) has base \(mu+x\).
The [linked Gauss chord estimate](MULTIVARIATE_LINKED_GAUSS_MATCHING_SEPARATION.md#angular-stability-near-a-mirrored-profile)
bounds the total log difference, simultaneously for
\(m<J\), by
\[
 \frac18\sum_{\rm matched}
    \frac{(J-1)|u-v|+|x-y|}{\min(u,v)^2}
              +\frac{J-1}{8}\sum_{\rm unmatched}u^{-1}.
\]
The unmatched sum includes slots of either profile.
If this bound satisfies either stated tolerance, the same all-degree
criterion applies. In a common-slope match the error
has no order multiplier. This allows asymmetric
profiles with no exact mirror rows and with arbitrary
reference curvature.

## Opposite-slot midpoint matching

When the numbers of positive and negative shared Gauss slots
are equal, pair all of them across orientations. Pairing can
use different slopes. For a pair \((u,x),(v,y)\) let
\(\Delta_m=m(u-v)+x-y\), \(u_0=\min(u,v)>0\), and use the
average base \(m(u+v)/2+(x+y)/2\) on both sides of the
reference. Such a reference is in the broader mirrored Gamma
class above, using a unit gap and positive real shape per slot.
It retains the two ordinary private units.

For a base \(A\) write \(f_k(A)=\log(A)_k-\frac{k}{m}\log(A)_m\).
For \(k\le m/2\) its antisymmetric derivative satisfies
\[
 0\le f'_k(A)-f'_{m-k}(A)
                 \le\frac{k(m-k)(m-2k)}{3A^3}.
\]
Indeed for \(k<m/2\) this difference is \(m-2k\) times
the difference between the full and middle uniform averages
of \(j\mapsto1/(A+j)\), over \(j=0,\ldots,m-1\)
and \(j=k,\ldots,m-k-1\), respectively.
The full uniform law is greater in convex order: the
intervals have the same center, and the full law's distance
from the center is stochastically larger. Subtracting
\(j^2/A^3\) from \(1/(A+j)\) gives a concave function.
The variance difference \(k(m-k)/3\) therefore yields
the upper bound. At \(k=m/2\) it is immediate.
Also \(|f''_k(A)|\le k(m-k)/A^3\), by the two-block
sum for the chord derivative and the Lipschitz bound
\(2/A^3\) on \(j\mapsto1/(A+j)^2\).

Taylor expansion at the average base and reflection from
\(k\) to \(m-k\) now give the pair's absolute normalized
log coefficient error at most
\[
 \frac{k(m-k)}{m^3 u_0^3}
           \left(\frac{|m-2k|\,|\Delta_m|}{6}
                              +\frac{\Delta_m^2}{4}\right).
\]
No slotwise equality of offsets is needed. In particular
take \(D= (J-1)|u-v|+|x-y|\) for each pair. At \(m<J\),
on the root disks used for the quadratic angular reserve,
the absolute evaluation error is bounded by
\[
 \frac12\sum_{\rm pairs} u_0^{-3}
                         \left(JD/6+D^2/4\right).
\]
This uses \(1\)-Lipschitzness of the exponential on
nonpositive log coefficients, \(|z|^m<3\), and
\(\sum_{k=0}^m k(m-k)=m(m^2-1)/6\).
Thus if this last bound is strictly below
\(1/(32\cdot10^9 J^2)\), Rouché and the quadratic
reference reserve prove all-degree angular separation
of the original ordinary profile.
For common-slope pairs with unequal offsets of bounded
size this requires slopes on the scale of \(J\), even
when there is no exact mirror and the complete-row counts
on the two sides differ. Exact mirror pairs can be
included in arbitrary number since their error is zero.
