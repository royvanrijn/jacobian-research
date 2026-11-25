# Effective subcritical thermal fold comparison

Use the fully shared chord-centered phase and the real center
\(q_*,w_*,a_2,g\) of the
[effective shared quartic window](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md).
On the subcritical side put
\[
 \delta=\max\rho-1<0,\qquad
 r=\sqrt{a_2/(2g)},\qquad c=-4g r^4,\qquad M=nc,\qquad
 b_j=r(w_*-w_j)/c,\quad j=n,N=n+1 .
\]
The limiting upper thermal fold has field \(b_c=2i/(3\sqrt3)\).

**Theorem.** On every prescribed positive rational band
\(W_-\le W\le W_+\) one can compute positive rational constants
\(\delta_0,r_0,b_0,M_0\) such that consecutive common roots are
impossible under
\[
  -\delta_0\le\delta<0,\qquad r\le r_0,\qquad
      M\ge M_0,\qquad |b_j-b_c|\le b_0\quad(j=n,N).
                                                               \tag{1}
\]
The constants are unevaluated. This excludes the full fixed
disk in the physical thermal field, including arbitrarily large
Airy coordinates and arbitrary integer degree difference.
Conjugation gives the lower fold with the same constants.
One may take \(b_0<|\Im b_c|/2\). If both fields are in the
union of the two disks, their equal imaginary parts put them
in the same disk.
Since \(c\asymp\delta^2\) with computable band constants,
the mass threshold can also be imposed by a computable lower
bound on \(n\delta^2\).

In the dimensionless convention
\(\ell=(-4g)^{-1/4}\), \(\ell_0=\ell^{4/3}\),
\(r_{\rm dim}=\sqrt{2|a_2|\ell_0}\), the physical length
in (1) is \(r=\ell_0r_{\rm dim}\).
Thus \(c=\ell_0r_{\rm dim}^4\), \(b_j=(w_*-w_j)/r_{\rm dim}^3\),
\(M_j=T_j^2\), and \(H_j=b_j M_j^{3/4}\).
An upper bound on \(\ell_0\) from the curvature band turns the
physical length cap \(r_0\) into a uniform dimensionless cap.

## The exact lattice error and a common scaled contour

The real-center and Taylor bounds from the effective quartic
window give \(0<g_-\le-g\le g_+\) and a fixed physical
complex disk about \(q_*\). They give
\(|q_*-q_m|\le C|\delta|\) for the maximizer \(q_m\) of \(\rho\).
Since \(\rho'(q_m)=0\) and \(\rho''\) is uniformly bounded,
\(\rho(q_*)-1=\delta+O(\delta^2)\). Thus \(a_2<0\),
\(r^2\asymp|\delta|\), and \(c\asymp r^4\), with computable
bounds as stated. Write \(\Phi_w(q)=F(q)-wq\), where \(F\)
is the chord-normalized phase with real entropy included.
The local exact phase is
\[
 \Phi_{w_j}(q_*+h)-\Phi_{w_j}(q_*)
          =(w_*-w_j)h+a_2 h^2+g h^4+h^5 E(h),
\]
with analytic bounds on \(E\) and its derivatives.

Choose a small fixed rational physical radius \(L\) so
\(|E(h)h|\le g_-/4\) for \(|h|\le2L\), and similarly choose
\(L\) so the derivative of \(h^5E(h)\) on the real interval
is at most \(g_-|h|^3\). The quartic term dominates these
errors. For \(b_j\) bounded, \(|w_*-w_j|=O(r^3)\).
Hence at \(h=\pm L/2\) the real phase has a fixed strict
loss and outward-decreasing derivative. The quadratic
term has the favorable sign on both sides. Since
\(F''=(\rho-1)/(q(1-q))<0\), the gap persists to the
lattice endpoints \(0,1\). All these gap constants
follow from strict inequalities using the fixed \(L,g_-\)
and a small enough \(r_0\).

At an integer order \(j\) the exact summand divided by the
constant factorial coefficient is
\(j^{-1/2}A_j(q)e^{j\Phi_{w_j}(q)}\) on the retained analytic
neighborhood, \(q=k/j\). Its interpolation is the Stirling
form of the exact Gamma quotient. On that neighborhood
\(A_j=A_0+O(1/|j|)\), analytic in \(q\) and in \(j\) on a fixed
narrow sector about the positive ray. The uniform band
bounds control all needed derivatives and give a positive
lower bound for \(|A_0|\) near \(q_*\). Indeed the row
arguments divided by their respective \(m_i\) stay in a
fixed right half-plane sector. The half-logarithms are
bounded in sum by a constant times the row curvature, and
the Stirling remainders are summed using \(\sum1/m_i\le W_+\).
The leading amplitude includes the bounded entropy factor
\((2\pi q(1-q))^{-1/2}\).
[DLMF Stirling remainder bounds](https://dlmf.nist.gov/5.11.ii)

For clarity, the cotangent reduction here has a fixed, effective
exponential error in the *unscaled* order. Retain the interval
\([q_*-L,q_*+L]\), adjusting initially to neighboring half-grid
endpoints at distance at most \(1/j\).
The entropy upper bound for the binomial coefficient and the
shared Stirling ratio bounds apply on the full exact grid,
so the discarded terms have the fixed outer phase loss.
On the retained interval apply the
[finite cotangent summation](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour)
at fixed heights \(\pm\nu\).
At the right and left walls the real loss persists for small
fixed \(\nu\) and small \(r\): the real phase increase is at most
\(C\nu^2+|\Im w_j|\nu\), by real analyticity.
On the horizontal sides the real phase relative to the center
is at most \(C r^3 L+C\nu^2+|\Im w_j|\nu\), by the negative
quadratic and quartic terms on the real diameter.
Choose \(\nu\) first so \(C\nu^2\) is less than both a fixed
fraction of the wall loss and \(\pi\nu/4\).
Then shrink \(r_0\) so the other terms obey the same strict
margins. The cotangent on half-grid walls has modulus at
most one and on the horizontal panels its nonconstant part
is \(O(e^{-2\pi j\nu})\). The constant parts give the real
integral. Moving the endpoints back to the fixed physical
interval costs only its strict exponential loss.

Now set \(q=q_*+rz\). The scaled phase is
\[
 \Psi_r(z,b)=-z^4/4-z^2/2+bz+r z^5E_r(z),
 \qquad E_r(z)=E(rz)/(-4g).
\]
Shift the finite interval to \(\Im z=1/\sqrt3\).
The two physical end connectors have the same fixed loss,
by further reducing \(r_0\) to bound their unscaled height
\(r/\sqrt3\). Let \(J(\mathcal M,b)\) be the entire retained
horizontal integral of \(A_{\mathcal M/c}(q_*+rz)
e^{\mathcal M\Psi_r(z,b)}\,dz\), where \(\mathcal M\) denotes
a complex interpolation of the mass.
At integer \(j\), writing \(M_j=jc\), the normalized exact
moment is this integral with an error
\[
                      O(r^{-1}j^C e^{-\gamma j}).           \tag{2}
\]
The normalization here divides the original polynomial by its
constant coefficient, by \(\sqrt j\), and by
\(r e^{j\Phi_{w_j}(q_*)}\).
Since \(M\ge1,\ c\le4g_+r^4\), one has
\(r^{-1}\le(4g_+n)^{1/4}\).
Thus (2) is \(O(e^{-\gamma'j})\) with computable
\(\gamma'>0\) after increasing a fixed terminal threshold.

## Constructive Airy chart and relative remainder

At \(z_c=i/\sqrt3,\ b_c=2i/(3\sqrt3)\) the model horizontal
real-phase loss is exactly \(-(\Re z)^4/4\) relative
to its value at \(z_c\).
In a small fixed central interval its cubic coordinate has
the expansion, with \(h=z-z_c\),
\[
   u=i3^{1/6}h(1-i\sqrt3 h/4)^{1/3}
       =i3^{1/6}h+3^{1/6}\sqrt3\,h^2/12+O(h^3).
\]
This follows from the stationary increment
\(-ih^3/\sqrt3-h^4/4\). For two sufficiently small fixed
real endpoints \(h=\pm h_0\), the images lie strictly in
the open right Airy decay sectors. Choose such a rational
\(h_0\) by the convergent binomial series or its explicit
cubic height at those endpoints.
In each sector connect that endpoint to a point of equal
radius on the central ray, then to a common positive
cutoff radius inside the chart. The real cubic phase on
these circular and radial joins is strictly negative.

The actual fold solves \(\partial_z\Psi_r=\partial_z^2\Psi_r=0\).
The model Jacobian in the two variables \(z,b\) is
invertible because \(\Psi_{zzz}=-2i\sqrt3\) and
\(\Psi_{zb}=1\) there. A fixed Newton contraction using
the model inverse Jacobian gives the actual fold and
its \(O(r)\) displacement on a fixed small polydisk.
The contraction bounds are obtained from the explicit
\(O(r)\) perturbation on larger scaled disks and Cauchy.
At the actual fold absorb the nonzero cubic unit into
a local coordinate \(v\), by the cube root near its
nonzero central value. Rouché gives a bounded inverse
on a fixed smaller disk.

The exact normal form
\[
                C_r(b)+u^3/3-\lambda_r(b)u                 \tag{3}
\]
then follows by the effective
[weighted Taylor norm contraction](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#a-uniform-cubic-chart-with-effective-constants).
In detail, for \(u=v+h\) the linear inverse removes
the constant and linear coefficients and divides the
remainder by \(v^2\). On the norm
\(|a|+R|\lambda|+R^2\|h\|_R\) it has norm one.
The nonlinear part \(vh^2+h^3/3-\lambda h\) has
Lipschitz constant at most \(4\eta/R^3+\eta^2/R^6\)
in a ball of radius \(\eta\), so one may take
\(\eta=R^3/100\). Bounds on the inverse coordinate on
a larger disk and the field derivative make the
analytic defect smaller than \(\eta/2\) on a fixed
field disk. This supplies uniform analytic coordinate
bounds and, by the linear term at the fold,
\(\lambda'_0(b_c)=i3^{-1/6}\).
On a smaller convex field disk the imaginary part of
\(\lambda_r'\) is bounded below, and its inverse has
imaginary derivative bounded above by a negative
constant. The inverse exists on one fixed disk about
zero by linear comparison on the perimeter.
All the smallness requirements follow from explicit
derivative bounds and positive margins.

The two strict Airy exit joins persist under these
perturbations and give the elementary Airy contour.
Outside the central interval on the truncated scaled
line, every bounded complement retains fixed real loss
relative to \(C_r(b)\), by the model quartic loss.
For the unbounded part up to \(L/r\), use large fixed
\(|x|\ge U\) on \(z=x+i/\sqrt3\).
The leading real term is \(-x^4/4+O(x^2+|x|+1)\);
the remainder \(r z^5E_r\) has modulus at most
\(C L |x|^4\) there. Decrease \(L\) at the start so
this is less than \(|x|^4/32\).
This gives a strict quartic tail loss, also relative
to the bounded center \(C_r(b)\).
The imaginary part on this tail is at most \(C|x|^4\).
Choose a fixed small relative complex mass disk around
each real \(M\). On the bounded contour pieces its
rotation preserves the strict real gaps by the bounded
imaginary parts; on the tail it preserves quartic loss.
It also keeps \(j=\mathcal M/c\) in the Stirling sector.
These tests give one computable disk fraction
\(\eta>0\) valid for all \(|\mathcal M-M|\le\eta M\).

Remove \(e^{\mathcal M C_r(b)}\mathcal M^{-1/3}\)
and the elementary contour normalization from \(J\).
Put \(z=\mathcal M^{2/3}\lambda_r(b)\), using the branch
positive for positive real \(\mathcal M\). The resulting
interpolation is
\[
 I(\mathcal M,b)=P(\mathcal M,b)\operatorname{Ai}(z)
   +\mathcal M^{-1/3}Q(\mathcal M,b)\operatorname{Ai}'(z)
   +E(\mathcal M,b),\qquad
 |E|\le C M^{-1}J_{\!A}(z),                                 \tag{4}
\]
where \(J_{\!A}=|\operatorname{Ai}|+
(1+|z|)^{-1/2}|\operatorname{Ai}'|\).
Here \(P\) is bounded above and away from zero, and
\(Q\) is bounded, with uniform analytic bounds and field
derivatives on smaller disks. Indeed divide the local
amplitude \(A_{\mathcal M/c}(q_*+r z(u,b))z_u(u,b)\)
by \(u^2-\lambda_r(b)\) using Cauchy interpolation.
The constant coefficient is uniformly nonzero at the
center and stays so after shrinking. The divided
remainder and its derivative are bounded by Cauchy
on fixed circles. Integrate by parts once, then apply
the [effective elementary Airy joint norm](AIRY_EFFECTIVE_ELEMENTARY_JOINT_NORM.md)
on the finite rays after scaling \(u\) by
\(\mathcal M^{1/3}\). The angle changes by a fixed small
rotation, and the lemma's size condition is enforced
by a fixed small bound on \(|\lambda_r(b)|\).
The strict-loss boundary terms, joins and exterior
pieces are absorbed by the explicit global lower bound
in that lemma, after shrinking this bound further.
All constants in (4) are computable from the fixed chart.

The exact integer-order moment and (4) differ relatively
to \(J_{\!A}\) by \(O(e^{-\gamma'' j})\).
In fact the normalization of (2) costs a fixed power of
\(M\) and \(\exp(O(M))\), the latter also bounding the
inverse joint norm for small \(\lambda_r\).
Require \(c=M/n\) small enough, by decreasing \(r_0\)
after the losses have been fixed. Thus those factors
are absorbed into \(\exp(-\gamma'j)\).
The complex-order interpolation has kept the whole
finite horizontal contour. Only the integer lattice
error and the physical endpoint connectors are omitted.

## Airy localization and the linked order increment

At a true root, (4) and the exponential lattice error
give
\[
    |\operatorname{Ai}(z)|\le C_* M^{-1/3}
                                      |\operatorname{Ai}'(z)|
\]
with computable \(C_*\). Fix this constant on the initial
chart and then further decrease the physical field disk.
For any \(p>0\), if
\(|\operatorname{Ai}(z)/\operatorname{Ai}'(z)|\le p\)
and \(p^2(|z|+2p)\le1/32\), there is an Airy zero
within \(2p\) of \(z\). In fact, the normalized solution
\(Y(v)=\operatorname{Ai}(z+v)/\operatorname{Ai}'(z)\)
of
\[
 Y(v)=Y(0)+v+\int_0^v(v-s)(z+s)Y(s)\,ds
\]
has norm at most \(4p\) on that disk, and its difference
from the affine part is at most \(p/4\). Rouché applies.
Use \(p=C_* M^{-1/3}\) and the small bound on \(|\lambda_r|\).
The Airy zeros are negative real and simple.
[DLMF Airy zero law](https://dlmf.nist.gov/9.9.i)

Thus \(z=a_k+O(M^{-1/3})\) and
\(|\Im\lambda_r(b)|\le C/M\) at each root.
The exact relation
\[
       b_N-b_n=(r/c)\Delta\log(N/n)\in\mathbb R
\]
and the effective lower bound for \(\Im\lambda_r'\) on
the common convex disk give \(|b_N-b_n|\le C/M\).
Therefore their two scaled Airy fields differ by at
most \(C_1 M^{-1/3}\), since \(M_N=M+c\) and \(c\le1\).
They have the same Airy label. To check this quantitatively,
at \(a_k\) use the same ODE integral equation with initial
data \(Y(0)=0,Y'(0)=1\). It gives
\(Y(v)=v[1+O((1+|a_k|)|v|^2)]\) and
\(Y'(v)=1+O((1+|a_k|)|v|^2)\).
On a disk of radius \(C_2 M^{-1/3}\) large enough
to include the preceding errors, the product
\((1+|a_k|)|v|^2\) is small after reducing the bound on
\(|\lambda_r|\), since
\(|a_k|M^{-2/3}\le\sup|\lambda_r|+O(1/M)\).
That disk contains no other Airy zero.

For this fixed label, make the substitution
\(b=\lambda_r^{-1}(\mathcal M^{-2/3}z)\) for
\(|\mathcal M-M|\le\eta M\) and
\(|z-a_k|\le R M^{-1/3}\), with a large fixed \(R\).
It stays in the chart after shrinking the allowed
physical field disk. Divide (4) by
\(P(\mathcal M,b)\operatorname{Ai}'(a_k)\).
The resulting analytic function is
\[
 Y(z-a_k)+\mathcal M^{-1/3}(Q/P)Y'(z-a_k)+O(1/M).
\]
The ODE bounds, with \(R\) chosen before the final decrease
of the \(\lambda_r\)-disk, make the first term dominate
on the boundary of this disk. Rouché gives a unique simple
root \(z=a_k+\chi_k(\mathcal M)\) with
\(|\chi_k|\le R M^{-1/3}\).
Cauchy's estimate on the smaller mass disk then gives
\(|\chi_k'|\le C M^{-4/3}\) for \(M\le\mathcal M\le M+c\).

For the comparison with actual roots, the normalized
\(z\)-derivative is within \(1/4\) of one on a smaller
disk containing both roots.
For the middle term above use
\((Q/P)_z=O(M^{-2/3})\) and
\(M^{-1/3}|Y''|=O((1+|a_k|)M^{-2/3})\);
for the last term use Cauchy on the larger disk.
An intermediate Rouché test keeps the analytic root
inside this smaller disk, and the actual root is there
by its localization.
The integer-order error divided by
\(|P\operatorname{Ai}'(a_k)|\) remains exponentially
small in \(n\) by the explicit joint lower bound at
\(a_k\) and \(M\le c n\) with \(c\) small.
Integration of the normalized derivative along the
segment to the analytic root therefore gives distance
\(O(e^{-\gamma'''n})\) in \(z\).
The analytic inverse \(\lambda_r^{-1}\) gives the same
bound for the distance in \(b\), up to a bounded factor.

Finally the analytic branch in the thermal field satisfies
\[
 \frac{d b_k}{d\mathcal M}
   =(\lambda_r^{-1})'
        \left(-\tfrac23 a_k\,\mathcal M^{-5/3}+O(M^{-2})\right).
\]
Its imaginary derivative has constant sign and modulus
at least \(c_1 |a_k| M^{-5/3}\) for \(M_0\) large:
the inverse derivative has strictly negative imaginary part,
and the negative Airy zeros are bounded away from zero
by the value at zero and the local analytic series.
Integrating over the actual increment \(c=M/n\) gives
an imaginary displacement at least
\(c_1 |a_k|\,n^{-1} M^{-2/3}\ge c_2' n^{-5/3}\),
using \(M\le n\).
The actual integer-order errors cannot cancel this
displacement. This contradicts the real exact field drift
without any bound or case split on \(\Delta\).

The choices occur in the following order. First fix the
physical interval and cotangent height from the band
derivative bounds, then the scaled contour, fold disks,
mass sector and fixed positive losses. Obtain the amplitude
and joint Airy constants on these fixed domains.
Choose the root-localization constant, the label and Rouché
radii, and then shrink the \(\lambda_r\)-disk to meet their
ODE conditions. Choose \(b_0,r_0\) so both physical fields
stay in that disk and so the lattice normalization leaves
a positive exponential loss. Finally increase \(M_0\) to
meet the finitely many bounds on the errors and the linked
increment. These operations give the constants of (1).
