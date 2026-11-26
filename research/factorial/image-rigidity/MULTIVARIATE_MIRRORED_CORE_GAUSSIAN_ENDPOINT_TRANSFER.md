# Angular separation for exact mirrored cores with Gaussian and diffuse tails

Let \(C_m(z)\) be the endpoint-normalized order-\(m\) factorial
polynomial of an ordinary core \(C\), obtained by a scalar and
positive change of variable. For \(\tau\ge0\) form its
physical Gaussian deformation
\[
 F_m(z)=\sum_{k=0}^m [z^k]C_m(z)\,
                        e^{-\tau k(m-k)/(2m)}z^k .
\]
The same \(\tau\) is used at both orders, so the angular heat
times are \(\tau/m\).

Suppose first that the core is mirrored. Write its shared pairs
as \((\ell_i,\ell_i+d_i),(\ell_i+d_i,\ell_i)\) with
\(\ell_i,d_i\ge1\), and write \(S_C\) for the sum over all unequal
rows \(d_i^2/\max(a_i,b_i)\). If the core is fully shared, put
\[
 V_C=2\sum_i d_i^2/\ell_i,\qquad W=1+V_C+\tau,\qquad
 Q(w)=e^{-\tau(1-w)/2}
           \prod_i \left(\frac{\ell_i+d_i w}{\ell_i+d_i}\right)^{d_i},
 \qquad q=Q(0).
\]
The following conditions suffice for strict interlacing of
\(F_n,F_{n+1}\) on the unit circle, in the angular ordering cut
at \(1\), at every \(n\ge1\):
\[
        Q(1/30001)\le10^{-6},\qquad Q(1/2)\le1/55,\qquad
                                      Wq\le10^{-4}.             \tag{1}
\]
Alternatively, a mirrored core containing a private pair has the
same conclusion if \(S_C+\tau\ge14\).
In both cases the circular cross-gap in root arguments at
\(n,N=n+1\) is greater than \(\pi/(20N^2)\).

In particular, a fully shared mirrored core with one pair
\(d_i\ge7\ell_i\) satisfies (1) for every \(\tau\ge0\), with
arbitrary other shared mirrored pairs. Every fully shared
mirrored core satisfies (1) for \(\tau\ge36\), including the
empty core.

## Endpoint comparison with the exact Gaussian clock

Write \(c_{m,k}=[z^k]F_m\), \(S=S_C+\tau\), and
\(V=V_C+\tau\) in the fully shared case.
The ordinary
[endpoint coefficient estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
gives the same bound with augmented \(S\), since the Gaussian
factor supplies precisely the additional exponential of curvature
\(\tau\). Similarly the exact row-product estimate for the
[hard mirrored core](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#a-criterion-using-all-the-endpoint-factors)
becomes, for \(1\le k\le m/2\), \(w=k/m\),
\[
 c_{m,k}^{1/k}\le \frac{eQ(w)}w,\qquad
 c_{m,k}\le \frac{(mq)^k}{k!}e^{V k^2/(2m)}.                  \tag{2}
\]
The first comes from multiplying the factorial ratios and the
exact Gaussian coefficient. For the second, use
\(\log(1+(d_i/\ell_i)w)\le(d_i/\ell_i)w\).
The augmented bounds \(S\le V\le S(S+2)\) hold because they hold
at \(\tau=0\), and \(S_C(S_C+2)+\tau\le(S_C+\tau)(S_C+\tau+2)\).

For the endpoint maximum in (2), set \(r_i=d_i/\ell_i\).
The derivative of \(\log(eQ(w)/w)\) multiplied by \(w\) is
\[
             -1+\sum_i\frac{d_i r_i w}{1+r_iw}+\frac{\tau w}{2},
\]
which increases in \(w\). Thus the maximum on a closed subinterval
is at an endpoint. Similarly \(m\) times the derivative of
\(\log(mQ(1/m))\) equals
\(1-\sum_i d_i r_i/(m+r_i)-\tau/(2m)\), which increases in \(m\).

Here are the checks that transfer the earlier analytic estimates
to the Gaussian deformation. Write \(N=n+1\) and use the symmetric
half \(H_m(z)=\sum_{k<m/2} c_{m,k}z^k+
\tfrac12\boldsymbol 1_{m\ {\rm even}}c_{m,m/2}z^{m/2}\).

- If \(N\le30001\), the two endpoint maxima above and (1) give
  \(c_{m,k}\le10^{-k}\) for \(1\le k\le m/2\), \(m\ge2\).
  If \(N>30001\) and \(S\ge2\log N+5\), the large-curvature
  geometric bound from the first estimate gives \(c_{m,k}\le8^{-k}\).
- Otherwise \(N>30001,\ S<2\log N+5\). We have \(V/(2n)<1/32\):
  with \(R=\log N>10\), the prior bound
  \(((2R+5)^2+(2R+5)/2)/(e^R-1)<1/32\) applies to \(V/(2n)\).
  Set \(D=Nq\). When \(D\le1/10\), (2) gives
  \(c_{m,k}\le(1/8)^k/k!\) for \(k\le7\), and testing
  \(k=8\) and \(k/m=1/2\) in its first bound gives
  \(c_{m,k}\le10^{-k}\) for \(8\le k\le m/2\).

In these regions \(H_m\) is zero-free on the unit disk and
\(\|zH_m'/H_m\|_{|z|\le1}<1/5\). For the geometric cap \(1/8\),
the ratio is at most \(4/21\), by summing the two geometric series.
The cap \(1/10\) is stronger. In the short \(D\) case, both
\(\|H_m-1\|\) and \(\|zH_m'\|\) are below \(3/20\) by summing
the factorial series and geometric tail. Their ratio gives
\(\|zH_m'/H_m\|<3/17\). Also \(H_1=1\).

For the long \(D>1/10\) region take \(K=\lceil32D\rceil+14\).
The inequalities \(q\le10^{-6}\), \(Wq\le10^{-4}\) imply
\(n\ge100000\),
\[
 D/n<1/720000,\quad WD/n<1/9000,\quad
 K\le182D,\quad WK/n<182/9000,\quad K+1<n/2 .
\]
The endpoint bounds at \(k=K\) and \(k/m=1/2\) in (2) give
\(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\); at the first endpoint
use \(eD e^{VK/(2n)}/K<1/8\).
For real \(m\in[n,N]\), the original factorial first coefficient
\(B_m^C\) is replaced by
\[
 B_m=B_m^C e^{-\tau(1-1/m)/2}.
\]
Its logarithmic derivative satisfies
\(0\le1/m-(\log B_m)'\le V/m^2\) by the original row-product derivative
bound and the added term \(-\tau/(2m^2)\).
Thus \(|(\log B_m)'|\le2/n\) and \(B_m/m\) is decreasing.
Equation (2) at \(m=n,k=1\) gives \(B_m\le D e^{V/(2n)}<101D/100\),
since \(V/n<1/900\).

In the
[endpoint interpolation](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation),
the logarithm of \(c_{m,k}k!/B_m^k\) gains exactly
\(\tau k(k-1)/(2m)\). This term contributes only to the second
finite difference, by \(\tau/m\), and its order derivative by
\(-\tau/m^2\). Both are absorbed by increasing \(V_C\) to \(V\)
in the derivative bounds (5)–(7) of that interpolation.
The scale \(\alpha=2\max(1,\max_i d_i/\ell_i)/n\) is unchanged
and \(\beta=8eW/n\) uses the augmented \(W\). The long \(D\)
inequalities give the same estimates
\(\alpha K<1/20,\ \beta K<1/2,\ \beta B_m<1/300\)
as in the
[paired long endpoint comparison](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#a-paired-estimate-for-the-two-long-endpoint-sums).
The connected-sum and two Cauchy tail estimates there thus apply
to these exact coefficients and their order derivatives. In
particular \(H_m\) has a logarithm zero at zero with
\[
 \|\log H_m\|_{|z|\le6/5}\le3(D+1),\qquad
 \|\log H_N-\log H_n\|_{|z|\le6/5}<1/55 .
\]
Cauchy's estimate gives on the unit circle logarithmic angular
derivatives at most \(15(D+1)<m/2\) for each half and less than
\(1/11\) for their difference. All of these bounds use the linked
clock in the first coefficient and the second edge weight.

For a mirrored core containing a private pair, the
[private-start factor](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
gives the bound \(e^{-(S_C-5)k/4}\) for \(k\le m/2\).
Multiplication by the Gaussian contributes at most
\(e^{-\tau k/4}\) in this range. For \(S_C+\tau\ge14\) the
result is below \(8^{-k}\), so the preceding geometric argument
applies without a cluster estimate.

## Strict interlacing and the angular reserve

On the circle, the centered value of \(F_m\) is
\(2|H_m(e^{i\theta})|\cos\psi_m(\theta)\), where
\(\psi_m(\theta)=m\theta/2-\phi_m(\theta)\) and
\(\phi_m=\Im\log H_m(e^{i\theta})\).
The logarithm is real on the real diameter, so \(\phi_m\) vanishes
at \(0,\pi\). The bounds above give \(0<\psi_m'<m\) and
\[
        |\phi_N(\theta)-\phi_n(\theta)|<2\theta/5
                          \qquad(0<\theta\le\pi).
\]
At a lower root with positive argument \(\theta\le\pi\),
\(\psi_n(\theta)\) is a root level \((k+\tfrac12)\pi\), and
\(\psi_N(\theta)-\psi_n(\theta)\) lies strictly between
\(\theta/10\) and \(9\theta/10\). Thus the upper phase at this
point is more than \(\theta/10\) from either neighboring root
level. This gives strict interlacing by phase monotonicity,
reflection and \(\psi_m(\pi)=m\pi/2\).
It also gives distance greater than \(\theta/(10N)\)
to any upper root argument in \([0,\pi]\).
The smallest positive lower root is at least \(\pi/(2n)\)
since \(\psi_n'<n\). Reflecting a root from the other half
to \([0,\pi]\) can only decrease its circular distance from a
fixed point in \([0,\pi]\). The claimed \(\pi/(20N^2)\)
cross-gap follows.

An exact fully shared mirrored core satisfying (1) at \(\tau=0\)
continues to satisfy it when \(\tau\) increases. The two \(Q\)
values decrease. Also \(q|_{\tau=0}\le10^{-6}\) implies
\(V_C\ge-2\log(q|_{\tau=0})>1\), so
\((1+V_C+\tau)e^{-\tau/2}\) decreases. The
[large-gap core estimate](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#a-criterion-using-all-the-endpoint-factors)
verifies (1) at zero for a pair with \(d_i\ge7\ell_i\).
For an arbitrary fully shared mirrored core,
\(Wq\le(S+1)^2e^{-S/2}\) since
\(\log(1+d_i/\ell_i)\ge d_i/(\ell_i+d_i)\).
For \(\tau\ge36\), this bound is at most \(37^2e^{-18}<10^{-4}\),
and the two \(Q\) evaluations meet (1) by the Gaussian factor alone.

## Ordinary asymmetric diffuse tails

Append arbitrary ordinary unequal rows to a mirrored core.
Set \(\tau=\sum_{\rm tail}d_i^2/\max(a_i,b_i)\) and
\(\eta=\max_{\rm tail}|d_i|/\max(a_i,b_i)\), with the empty tail
allowed. If the exact Gaussian reference for this \(\tau\)
satisfies either criterion above, then the full ordinary profile
has disjoint consecutive root arguments at every order for
all \(\eta\le\eta_*\), for one computable absolute
\(\eta_*>0\). No sign or balance condition is imposed on the tail.

For this claim, take \(J\ge4\) at least the
[effective absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
handles the exact ordinary case of \(S_C+\tau\ge27\), so
assume \(\tau<27\). Let \(0<\eta_0\le1/2\) be a rational
tolerance from the
[absolute ordinary diffuse theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md).
If all core relative gaps are at most \(\eta_0\) and
\(\eta\le\eta_0\), that theorem handles the full profile.
Otherwise a core pair of relative gap greater than \(\eta_0\)
provides a Gegenbauer factor of parameter
\(1+m\ell/d\le1+m/\eta_0\).
The
[Gegenbauer mesh and circular composition bound](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles)
therefore bounds the chord gap in \(F_m\) below by
\(g=\min(1,4/[J(2+1/\eta_0)])\) for \(2\le m<J\).
Here the Gaussian factor is itself circle-rooted and composition
preserves the lower mesh bound of the indicated core factor.

Take \(r=\min(g/4,1/(100J^2))\) and
\(\eta_*=\min(\eta_0,r(g/2)^{J-1}/(27J\,3^J))\).
The
[diffuse tail coefficient bound](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md#a-comparison-at-exactly-the-same-signed-center)
for the endpoint-normalized full \(P_m\) gives
\(|P_m(z)-F_m(z)|\le(m\tau\eta/4)(1+|z|)^m\).
It uses the exact factorial core in both polynomials.
On every disk boundary of radius \(r\) about an \(F_m\) root,
this is strictly below \(r(g/2)^{J-1}\le |F_m(z)|\).
Rouché puts a root of \(P_m\) in each disk, hence every full root
argument moves by at most \(2r\).
For terminals \(3\le N<J\) the two errors total at most \(4r\),
smaller than the Gaussian reference gap. The first-to-second
pair separates by positivity of the ordinary quadratic moment
on the real line. Terminals \(N\ge J\) follow from the absolute
ordinary cutoff.

## Absolute curvature for a general ordinary core plus Gaussian

The same augmentation extends
[the absolute curvature 27 theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
to a possibly asymmetric ordinary core, with private rows
allowed. If \(S_C+\tau\ge27\), then the exact Gaussian-deformed
moments \(F_n,F_{n+1}\) have disjoint root arguments at every order.

To check the hypotheses of that proof, use \(S=S_C+\tau\)
and \(V=V_C+\tau\), where \(V_C\) there includes the private
terms. The bound \(S\le V\le S(S+2)\) still holds and
the geometric coefficient estimate gains exactly the
Gaussian curvature \(\tau\). In a private-start half the
coefficient estimate becomes \(e^{-(S-5)k/4}\) as above.
For a positive starting endpoint, the log interpolation
again gains \(\tau k(k-1)/(2m)\). Its finite-difference and
order-derivative contributions are absorbed by the augmented
\(V\). The absolute curvature proof's bound for the first
coefficient derivative gains at most \(\tau/(2m^2)\), which
is bounded by the increase from replacing \(V_C,S_C\) by
\(V,S\) in \((5V+S(3+\log N))/m^2\). The quantities \(\alpha\)
and \(\beta\) retain the required relations because the row
ratio maximum is unchanged and \(W=1+V\) increases.
All subsequent bounds on \(D,K,\Xi\) use only \(S,V\), the
first coefficient and its derivative, the coefficient tails,
and these finite-difference bounds. Thus the same angular
argument applies. For the first-to-second pair, the Gaussian
factor realizes positive lognormal averaging of the even
moment. That average is strictly positive on the real axis
for \(\tau>0\), and the exact ordinary pair separates at
\(\tau=0\) when \(S_C\ge27\).


For mirrored cores one also obtains strict interlacing and a cross-gap
greater than \(\pi/(240N^2)\) throughout \(S_C+\tau\ge27\).
Use the symmetric half with half of the middle term. In the two
geometric regions of the absolute curvature proof, its logarithmic
angular derivative on the unit circle has modulus less than
\((1/5)/(1-37/200)=40/163\). In the large-curvature geometric branch
\(S\ge2\log N+5,\ N\ge256\), use \(4/21\) as above.
In the cluster region the same symmetric half agrees through \(K\)
with the connected logarithm approximation, and its tail satisfies
the same bound. Thus estimate (3) of the absolute curvature proof
gives the individual and paired logarithm bounds for that half.
They give \(|\phi_m'|\le15(D+1)<m/2\) and paired derivative less than
\(1/11\), by Cauchy's estimate on radius \(6/5\) and the bound
\(D/n<1/720000\). Hence in all regions
\(|\phi_N-\phi_n|<(80/163)\theta\) and \(0<\psi_m'<m\).
The margin in the phase argument above is now
\((1/2-80/163)\theta=3\theta/326>\theta/120\),
giving the stated angle bound in the same way. The short
geometric estimates cover \(m=2\) and \(H_1=1\).
