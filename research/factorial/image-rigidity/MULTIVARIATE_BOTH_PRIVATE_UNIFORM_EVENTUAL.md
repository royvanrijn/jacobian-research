# A profile-uniform terminal cutoff when both monomials have private variables

Let \(\Lambda(X^\gamma)=\prod_i\gamma_i!\) be the ordinary
factorial functional. There is an absolute \(N_0\) such
that, whenever each of two distinct monomials has a
variable absent from the other, their consecutive moment
polynomials at orders \(n,n+1\) are coprime for \(n\ge N_0\).
The statement is uniform in the dimension and all exponents.

When both private degrees are at least two, the
[all-order private-degree theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two)
applies. When one private degree is one and the other is
at least two, the
[private-linear tail theorem](MULTIVARIATE_ONE_PRIVATE_LINEAR_TAIL.md)
gives a uniform cutoff. Here we treat the remaining case,
where each private degree is one. No numerical value of
the absolute cutoff is calculated.

## Normalization and coefficient estimates

Delete equal rows. When each private degree is one,
the two private factorials contribute \((n-k)!k!\);
after multiplying by \(\binom nk\) they cancel the
binomial exactly. Thus the endpoint-normalized polynomial
is
\[
 P_n(z)=\sum_{k=0}^n c(n,k)z^k,\qquad
 c(n,k)=\exp L_n(k/n),\qquad c(n,0)=c(n,n)=1 ,
 \tag{1}
\]
where \(L_n\) is the endpoint-subtracted Gamma log for
the shared rows alone. For the original parameter \(t\),
\(z=e^{\mu_n}t\) with
\[
 \mu_n=\frac1n\sum_{\rm shared}
       \log\frac{(n b_i)!}{(n a_i)!}.
\]
Within each consecutive pair the shared profile is fixed,
but it may depend on \(n\). Let
\[
 d_i=b_i-a_i,\quad m_i=(a_i+b_i)/2,\quad
 W=\sum_{\rm shared}d_i^2/m_i,\qquad
 \Delta=\sum_{\rm shared}d_i .
\]
If there are no varying shared rows, (1) is the geometric
sum and its consecutive pairs are coprime. Hence assume
\(W>0\).

A direct bound valid for all profiles is
\[
 c(n,k)\le \exp\!\left(-\frac{W}{8}\min(k,n-k)\right).
 \tag{2}
\]
Indeed the second coefficient difference of a shared row
is at least \(d_i^2/(n\max(a_i,b_i))\) on the log scale.
Subtract the corresponding quadratic with zero endpoint
values and compare with the endpoint chord, as in the
[fixed-profile proof](MULTIVARIATE_TWO_MONOMIAL_EVENTUAL_BORDER.md).
Then use
\(\sum d_i^2/\max(a_i,b_i)\ge W/2\) and
\(k(n-k)/n\ge\min(k,n-k)/2\).

For bounded \(W\), the
[varying-profile estimates](MULTIVARIATE_VARYING_PROFILE_PHASE.md)
give, for real order \(x\ge1\),
\[
 L_x=xG+H+R_x,\qquad
 \mu_{n+1}-\mu_n=\Delta\log(1+1/n)+O(W/n^2).
 \tag{3}
\]
The functions \(G,H,R_x\) vanish at the endpoints \(0,1\).
On \([0,1]\), \(G'',H''=O(W)\),
\(R_x''=O(W/x)\), and
\(\partial_x R_x''=O(W/x^2)\). The bounds for \(R_x\)
follow from the trigamma integral and its order
derivative in that note. The endpoint conditions give
the same bounds on the first derivatives and values,
including for \(\partial_x R_x\).
As \(W\to0\), the row ratios
\(|d_i|/(2m_i)\le W/2\) give
\[
 G''=W+O(W^2),\qquad
 G(u)=-\tfrac W2 u(1-u)+O(W^2) \quad(0\le u\le1).
 \tag{4}
\]

We record two real coefficient estimates used below.
For \(0\le k\le x\), use the smooth interpolation
\(c(x,k)=\exp L_x(k/x)\), and put
\(E(x,k)=\partial_x\log c(x,k)\), with \(k\) fixed.
Equations (3) give
\[
 |E(x,k)|\le C W(k^2+k)/x^2,\qquad
 |\partial_k E(x,k)|\le C W(k+1)/x^2.                     \tag{5}
\]
Indeed, for \(u=k/x\), the main term of \(E\) is
\(G(u)-uG'(u)=-\int_0^u vG''(v)\,dv\).
The \(H\) term is \(-u H'(u)/x\), and the terms
from \(R_x\) and its order derivative satisfy the
bounds above. Differentiating in \(k\) multiplies
the corresponding \(u\) derivative by \(1/x\).
Both estimates also apply to
\(c_{\rm R}(x,k)=\exp L_x(1-k/x)\), by reversing
the profile.

The second estimate concerns small \(W\) and
\(k\le 3x/4\). There is an absolute \(c_0>0\) with
\[
 c(x,k)\le e^{-c_0 W k}.
 \tag{6}
\]
To see this for real \(x\), use
\(\psi'(1+x y)\ge 1/(1+xy)\) for \(y>0\) in
\(L_x''\). Since \(y=a_i+d_i u\le2m_i\),
\(m_i\ge1\), one obtains \(L_x''\ge xW/3\).
Endpoint subtraction gives
\(L_x(u)\le -xW u(1-u)/6\).
The reversed estimate is identical.

## An angular comparison for split polynomials

Suppose an endpoint-normalized polynomial is written
\(P_m(z)=L_m^*(z)+z^m U_m^*(1/z)\), with both
half-polynomials positive on the positive real axis.
Suppose they are nonzero throughout the two disks
containing the normalized roots and their reciprocals.
For a common original root, the two normalized roots
have the same argument \(\theta\in(0,\pi]\), since
the two scalings are positive. Write \(r_m\) for their
respective radii and \(\varphi_m(\theta)\) for the
continuous argument of
\[
 L_m^*(r_m e^{i\theta})/U_m^*(r_m^{-1}e^{-i\theta})
\]
that is zero at \(\theta=0\). At a root
\(m\theta=(2j_m+1)\pi+\varphi_m(\theta)\).
Thus it suffices that
\[
 \sup_{[0,\pi]}|\varphi_{n+1}'-\varphi_n'|<1.              \tag{7}
\]
Indeed the two root equations imply
\(|\theta-2\pi(j_{n+1}-j_n)|<\theta\), impossible on the
indicated angular interval. The phase derivatives are
the real parts of the sums of the two logarithmic
derivatives \(z\,\partial_z\log L_m^*\),
\(w\,\partial_w\log U_m^*\) at their respective arguments.

## Curvature bounded away from zero

Suppose contrary to the theorem there is a sequence of
common roots at increasing \(n\). Pass to a subsequence
on which \(W\) tends to a finite nonzero value, tends to
infinity, or tends to zero. In this section we exclude
the first two cases.

In both cases split (1) at \(\lfloor m/2\rfloor\) and
reverse the upper half. If \(W\to\infty\), (2) gives
uniform convergence of both halves and their derivatives
to \(1\) and \(0\), respectively, on any fixed disk,
say \(|z|\le2\). On \(|z|=1/2\) the lower half dominates
\(z^mU_m^*(1/z)\) for large \(m\), and reversal gives
the matching outer bound \(|z|<2\) by Rouché's theorem.
Both logged phase slopes in (7) tend uniformly to zero.

If \(0<\inf W\le\sup W<\infty\), put
\(\rho_{\rm L}=e^{G'(0)}\), \(\rho_{\rm R}=e^{-G'(1)}\).
There is a fixed \(q_0<1\) bounding these parameters
and the rate in (2) above. In fact \(G''\) is comparable
to \(W\) uniformly in the real interval by bounded \(W\).
For each fixed coefficient index, (3) gives convergence
of the two halves' coefficients to the corresponding
powers of \(\rho_{\rm L},\rho_{\rm R}\), uniformly in
the profiles. Thus on \(|z|\le R\), for some fixed
\(1<R<1/q_0\), the halves and their derivatives differ
uniformly by \(o(1)\) from
\((1-\rho_{\rm L}z)^{-1}\) and
\((1-\rho_{\rm R}z)^{-1}\) and their derivatives.
Use (2) to sum the geometric majorants.
Both model functions have uniform positive lower
and finite upper bounds in modulus there.
Rouché and reversal give \(1/R<|z|<R\).
At a root \(r_m^m=|L_m^*|/|U_m^*|\), so \(r_m\to1\)
uniformly. The models use exactly the same
\(\rho_{\rm L},\rho_{\rm R}\) for both orders of the
pair, hence their angular log derivatives differ by
\(o(1)\) at the two normalized radii.
This proves (7) and excludes the second case too.

## Small curvature with \(nW\) unbounded

Now \(W\to0\) and, after a subsequence,
\(\mathcal L=nW\to\infty\).
Use a fixed smooth cutoff \(\chi\) with
\(0\le\chi\le1\), \(\chi(u)=1\) on \(u\le1/4\),
\(\chi(u)=0\) on \(u\ge3/4\), and
\(\chi(u)+\chi(1-u)=1\).
The half-polynomials now have coefficients
\(\chi(k/m)c(m,k)\) and \(\chi(k/m)c_{\rm R}(m,k)\).
The decomposition in (7) holds exactly.
The cutoff allows comparison of adjacent orders
without a discontinuous change at the midpoint.

Write \(\rho=\rho_{\rm L}=e^{G'(0)}\) and
\(\widetilde L_x(z)=(1-\rho z)L_x^*(z)\), allowing
real \(x\) in the interpolation.
The symmetric notation applies on the right.
For a small fixed \(b>0\), on \(|z|\le e^{bW}\)
one has, with \(D_z=z\partial_z\),
\[
 \begin{aligned}
 |D_z^j(\widetilde L_x-1)|
     &\le C/(x W^{j+1}) &&(0\le j\le2),\\
 |D_z^j\partial_x\widetilde L_x|
     &\le C/(x^2 W^{j+1}) &&(0\le j\le1),
 \end{aligned}                                               \tag{8}
\]
uniformly for \(x\in[n,n+1]\) and large \(\mathcal L\).
Here is a direct coefficient verification.
By (3) the increment
\(L_x(k/x)-L_x((k-1)/x)-G'(0)\) is \(O(Wk/x)\).
Consequently, before cutting off,
\[
 |c(x,k)-\rho c(x,k-1)|
       \le C e^{-c_0'W k} W(k+1)/x
\]
on the cutoff support, by (6).
Taking the first difference of the cutoff introduces
only \(O(e^{-c_0'W k}/x)\) on its transition region.
These are the coefficients of \(\widetilde L_x-1\).
Summing after multiplying by \(k^j e^{bW k}\)
gives the first line of (8), by
\(\sum_{k\ge1} k^j e^{-c_0'W k/2}=O(W^{-j-1})\).

For the second line, the order derivative of the
uncut difference is
\[
 (c_k-\rho c_{k-1})E(x,k)
        +\rho c_{k-1}(E(x,k)-E(x,k-1)).
\]
By (5) its modulus is bounded by
\(C e^{-c_0'W k}\{W^2(k+1)^3/x^3+W(k+1)/x^2\}\).
The terms produced by the cutoff and its derivatives
add only
\(C e^{-c_0'W k}\{W(k+1)^2/x^3+W(k+1)/x^2+x^{-2}\}\).
To see this, use first cutoff differences \(O(1/x)\),
order derivatives \(O(1/x)\) on their transition support,
and mixed index-order differences \(O(1/x^2)\).
Geometric summation gives the second line of (8),
using \(xW\ge1\). Extension by zero beyond the cutoff
support is smooth as \(x\) varies.

Use \(R=e^{bW}\), decreasing \(b\) if necessary.
Both denominators \(1-\rho_{\rm L} z\) and
\(1-\rho_{\rm R} z\) are nonzero on this disk since by
(4)
\[
 -\log\rho_{\rm L}=\tfrac W2+O(W^2),\qquad
 -\log\rho_{\rm R}=\tfrac W2+O(W^2).
\]
The tilded halves are nonzero there by (8).
On \(|z|=1/R\), the moduli of the model denominators
for \(z\) and \(1/z\) are comparable to
\(W+|1-e^{i\theta}|\), uniformly in \(\theta\).
Thus the ratio of the two half-polynomial moduli is
bounded, and \(R^{-n}\to0\) gives the inner bound by
Rouché. Reversal gives the outer bound.
At a root, comparing these model denominators with
\(\rho_{\rm L}-\rho_{\rm R}=O(W^2)\) gives
\[
 n|\log r_n|
 \le C\bigl(\mathcal L^{-1}+W+|\log r_n|/W\bigr).
\]
For this estimate each denominator has modulus at
least a fixed multiple of \(W\), and the two tilded
relative errors are \(O(\mathcal L^{-1})\).
Absorbing the last term gives
\(n|\log r_n|=o(1)\), also at \(n+1\).

For a common root, (3) and the radius estimates force
the integer \(\Delta\) to be zero eventually. Its
two normalized radii then differ on the log scale by
\(O(W/n^2)\). By (8), the change of the logged
tilded derivative at the same argument is
\(O((nW)^{-2})\), since the tilded halves
are bounded away from zero. The geometric denominator
is unchanged in this comparison. Moving the radius
by \(O(W/n^2)\) adds \(O(1/(n^2W))\), because the
second logarithmic derivative of each half is
\(O(W^{-2})\) on the comparison disks.
This follows from the geometric factor and the first
line of (8) through \(j=2\).
These estimates give (7) and exclude this regime.

## Small curvature with bounded \(nW\)

It remains to assume \(W\to0\) and \(\mathcal L=nW\)
bounded. Pass to a subsequence with
\(\mathcal L\to\Lambda\ge0\).
Here \(c(n,k)\), viewed as samples of a function of
\(k/n\), has uniform bounds through its second
derivative. This follows from (3).
Put
\[
 D_n(z)=\sum_{k=1}^n [c(n,k)-c(n,k-1)]z^k.
\]
Then
\[
 (1-z)P_n(z)=1-z^{n+1}+D_n(z).                            \tag{9}
\]
The first differences are \(O(1/n)\). For
\(|z|\le e^{-B/n}\) with a sufficiently large fixed
\(B\), their geometric sum and \(|z|^{n+1}\) are
together less than one. Equation (9) excludes roots
there. Reversal gives \(|\log r_n|\le B/n\) at all
roots.

For a root with \(n|1-z|\to\infty\), sum the differences
in \(D_n\) once more:
\((1-z)D_n\) has endpoint terms \(O(1/n)\)
and second differences \(O(1/n^2)\).
On the known annulus its modulus is \(O(1/n)\).
Thus \(D_n(z)=o(1)\) and
\(z^{n+1}=1+o(1)\), so \(n\log r_n\to0\).
If \(n|1-z|\) stays bounded, write \(z=e^{h/n}\) with
\(h\) bounded. The Riemann sum in (1) tends to
\[
 F_\Lambda(h)=\int_0^1 e^{-\Lambda u(1-u)/2+hu}\,du
       =e^{h/2-\Lambda/8}G_\Lambda(h),\qquad
 G_\Lambda(h)=\int_{-1/2}^{1/2} e^{\Lambda x^2/2+h x}\,dx .
 \tag{10}
\]
We next verify the model facts needed for this limit.

Every zero of \(G_\Lambda\) is imaginary and simple,
for every \(\Lambda\ge0\). For the location, sample the
weight \(e^{\Lambda x^2/2}\) on an odd-degree symmetric
grid to approximate the integral. The polynomial in
\(e^{h/n}\) so obtained is \(A(z)+z^n A(1/z)\)
with positive strictly decreasing coefficients in
the lower half \(A\), if \(\Lambda>0\).
Such an \(A\) has no zero in the closed unit disk:
in \((1-z)A(z)\) the positive constant dominates
the sum of the nonconstant moduli for \(|z|<1\).
On the circle equality would require \(z=1\),
using the strict first difference, and \(A(1)>0\).
The quotient \(z^n A(1/z)/A(z)\) has modulus one
on the circle and is holomorphic and nonconstant
inside, hence has modulus less than one inside.
Thus the full polynomial has no interior zero,
and reciprocity excludes exterior zeros.
The uniform Riemann convergence on compact \(h\)-sets
gives the location by Hurwitz. At \(\Lambda=0\)
the assertion is immediate from the hyperbolic sine.

For simplicity use
\[
 \Lambda G_\Lambda'(h)+hG_\Lambda(h)
      =2e^{\Lambda/8}\sinh(h/2),\qquad
 \partial_\Lambda G_\Lambda=\tfrac12 G_\Lambda'' .         \tag{11}
\]
The first identity is integration of a total derivative,
and the second differentiation under the integral.
A repeated zero at \(\Lambda>0\) must have
\(\sinh(h/2)=0\), and differentiating the first
identity shows it is exactly double.
For \(f(y,\Lambda)=G_\Lambda(iy)\) at such a zero
\((y_0,\Lambda_0)\), the scaled \(f(y_0+\epsilon z,\Lambda_0-\epsilon^2)\)
divided by \(\epsilon^2\) tends to
\(\tfrac12 f_{yy}(y_0,\Lambda_0)(z^2+1)\).
Rouché would give nonreal zeros at a positive
parameter slightly below \(\Lambda_0\), contradicting
the location just proved.

For \(\Lambda>0\) every positive zero \(y\) of
\(f(y,\Lambda)\) lies in a left half-period
\(((2j-1)\pi,2j\pi)\).
No zero can lie at \(2j\pi\) by (11) and simplicity.
At \((2j-1)\pi\) the cosine integral has the sign
of its last lobe, by grouping successive full
lobes of alternating signs with strictly increasing
absolute masses under the increasing weight.
The initial half-lobe has smaller mass than the
next full lobe. Thus these points are also excluded.
At \(\Lambda=0\), the positive zeros are \(y_j=2j\pi\),
and implicit differentiation of
\(f(y,0)=2\sin(y/2)/y\) using the heat identity gives
\(y_j'(0)=-1/y_j<0\).
A zero at positive parameter can be continued back
to zero parameter in its bounded half-period:
simplicity continues it at each positive parameter,
and the forbidden endpoints keep it bounded.
Its only possible limit is a zero at parameter zero.
The initial leftward motion then proves the
left half-period assertion.

It follows at every zero of \(G_\Lambda\) that
\[
 \left(h\partial_h+\Lambda\partial_\Lambda\right)
                         G_\Lambda(h)\ne0.               \tag{12}
\]
For \(\Lambda=0\) use simplicity and \(h\ne0\).
For positive \(\Lambda\), (11) shows that at a zero
\(h=iy\) the expression equals
\[
 \frac{e^{\Lambda/8}}2
       \left(\cos(y/2)-\frac{2y}{\Lambda}\sin(y/2)\right).
\]
For positive \(y\) the two trigonometric terms have
opposite signs in the left half-period; negative \(y\)
follows by evenness.

We return to the supposed common root. If the scaled
angles are unbounded, the high-angle estimate from
(9) gives \(n\log r_n\to0\) at both orders. If the
scaled angles are bounded, (10) and the imaginary
model zeros give the same conclusion by subsequential
limits of bounded \(h\). Hence \(\Delta=0\) eventually
by (3).

In the high-angle case the movement of the argument
between the two normalized polynomials is radial with
\(\log(z_{n+1}/z_n)=O(W/n^2)=O(n^{-3})\).
The mixed estimate (5) with \(W=O(1/n)\) shows that
for \(1\le k\le n\) the adjacent change of
\(c(n,k)-c(n,k-1)\) is \(O(n^{-2})\).
Indeed its order derivative is
\((c_k-c_{k-1})E(x,k)+c_{k-1}(E(x,k)-E(x,k-1))\).
The extra final difference at \(n+1\) is \(O(1/n)\).
These bounds and the radial motion give
\(D_{n+1}(z_{n+1})-D_n(z_n)=O(1/n)\) on the annulus.
Subtracting the two identities (9) at the roots
therefore gives
\((z_n-1)z_n^{n+1}=O(1/n)\), impossible when
\(n|1-z_n|\to\infty\).

Finally suppose the scaled angles are bounded. Put
\(h=n\log z_n\). The adjacent scaled parameter is
\(h'=h+h/n+O(n^{-2})\).
The trapezoidal estimate, with a uniform second
derivative bound, gives
\[
 P_n(e^{h/n})/n
   =\int_0^1 e^{L_n(u)+hu}\,du
          +\frac{1+e^h}{2n}+O(n^{-2}).
\]
In the difference at the two orders, the explicit
endpoint term changes by \(O(n^{-2})\).
By (3) the change of the integral, multiplied by \(n\),
tends along a subsequence with \(h\to h_0\) to
\[
 (h_0\partial_h+\Lambda\partial_\Lambda)F_\Lambda(h_0).
\]
Indeed its interpolated exponent changes by
\(G(u)+hu/n+O(n^{-2})\), with
\(nG(u)\to-\Lambda u(1-u)/2\).
At a common zero the expression would vanish, with
\(F_\Lambda(h_0)=0\) by (10).
Removing the nonzero scalar factor in (10) leaves
exactly the expression in (12), a contradiction.

Thus no sequence of common roots with unbounded
terminal order exists, including when the exponent
profile changes at every order. The characteristic-zero
claim follows from rational Bézout identities.
