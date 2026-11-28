# The preferred rotations of three circular slices can be incompatible

An angular pairing can be optimized separately for two prescribed
mass vectors. A third circle introduces a constraint because the
three relative rotations must add to zero. The individually optimal
rotations can violate this constraint even when all three mass vectors
are arbitrarily close to equal.

Fix a compact interval \(K\subset(0,\infty)\) of angular parameters.
There are three circular mass vectors of size \(O(\epsilon)\) from
\((1/3,1/3,1/3)\) whose average is \((1/3,1/3,1/3)\), such that every
ordinary realization of the three slices loses at least
\(c_K\epsilon^4\) in the sum of the three pairings compared with
optimizing the pairs separately. The three parameters may differ
within \(K\). In the planar Gaussian problem, placing these vectors
on three disjoint sets of equal radial measure in a compact annulus
gives a balanced mass profile with a strictly positive gap from the
integrated exact pair envelope.

## Preferred relative rotation

Write \(a_i=2\pi i/3\) for \(i=0,1,2\), with indices taken modulo
three and the endpoints lifted periodically. For \(x\in\mathbf
1^\perp\) small, let \(P(x)\) have endpoints \(a_i+x_i\) and let
label \(i\) be the arc ending at endpoint \(i\). The circular
masses are
\[
             p_i(x)=\frac13+\frac{x_i-x_{i-1}}{2\pi}.
\]
This map is invertible on the mean-zero increments. Let
\(F_\kappa(x,y,u)\) be the pairing of \(P(x)\) rotated by \(u\)
with \(P(y)\), against the normalized kernel
\(k_\kappa(t)=e^{\kappa\cos t}/I_0(\kappa)\) and angular
measure \(dt/(2\pi)\). Use \(F_\kappa(A,B)\) for the same pairing
on arbitrary labeled slices. The exact paired maximum for widths
\(p(x),p(y)\) is attained by consecutive arcs in a common cyclic
order for all sufficiently small \(x,y\), uniformly on \(K\).
[Exact pair envelope near balanced masses, JA2](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md#exact-pair-envelope-near-balanced-masses)

At \(x=y=0\) the relative rotation \(u=0\) is the unique
maximizer modulo \(2\pi\). In fact the pairing as a function of
\(u\) has nonnegative cosine coefficients, since every angular
Fourier eigenvalue of \(k_\kappa\) is positive. Its first
coefficient is strictly positive. Put
\(\Delta_\kappa=k_\kappa(0)-k_\kappa(2\pi/3)>0\). Differentiating
once at each moving endpoint gives
\[
      \partial_u^2 F_\kappa(0,0,0)
        =-\frac{6\Delta_\kappa}{(2\pi)^2}.
\]
Indeed the arc ending at an endpoint gains mass as the endpoint
moves, and the next arc loses the same mass. The derivative of
their score difference at a reference endpoint has value
\(-2\Delta_\kappa/(2\pi)\). The extra factor \(1/(2\pi)\)
comes from the moving mass. This agrees with the joint Hessian
in (JA3).

The implicit function theorem and compactness give a unique
optimal relative angle \(u_\kappa(x,y)\) near zero modulo \(2\pi\),
smooth in \(x,y\) and uniform on \(K\). Global uniqueness follows
by excluding a fixed neighborhood of zero using the strict
maximum at \(x=y=0\), then using strict angular concavity near
zero. Reversing both cyclic orders changes the sign of the
preferred relative angle.

The mixed endpoint derivative is
\[
  (2\pi)^2\partial_{x_i}\partial_{y_j} F_\kappa(x,y,u)
       =(d_i\cdot d_j)\,k_\kappa(a_i+x_i+u-a_j-y_j),
\]
where \(d_i=e_i-e_{i+1}\) describes the transfer at endpoint
\(i\). Thus the mixed quadratic term of
\((2\pi)^2\partial_u F_\kappa(x,y,0)\) is \(x^{\mathsf T}B_\kappa y\)
with
\[
       (B_\kappa)_{ii}=0,\qquad
       (B_\kappa)_{ij}=-k_\kappa'(a_i-a_j)\quad(i\ne j).
\]
The terms in \(x\) alone through second degree vanish. Indeed
\(\partial_u F_\kappa(x,0,0)\) is invariant under cyclic
permutation of endpoints and changes sign under reflection,
which acts on \(x\) by a signed permutation. The only
cyclically invariant linear form on \(\mathbf1^\perp\) is zero
and the only cyclically invariant quadratic forms there are
multiples of the squared Euclidean norm. The same argument
applies to \(y\) alone. Consequently
\[
    u_\kappa(x,y)=
         \frac{x^{\mathsf T}B_\kappa y}{6\Delta_\kappa}
         +O_K((|x|+|y|)^3).                                    \tag{RC1}
\]
Here the implicit equation first gives that \(u_\kappa=O_K((|x|+|y|)^2)\).
Substitution into its Taylor expansion gives (RC1) with the
stated remainder. Explicitly \(B_\kappa=b_\kappa M\), where
\[
  b_\kappa=\frac{\sqrt3}{2}
                \frac{\kappa e^{-\kappa/2}}{I_0(\kappa)}>0,\qquad
  M=\begin{pmatrix}0&-1&1\\ 1&0&-1\\ -1&1&0\end{pmatrix}.
\]

Take \(v=(1,-1,0)\), \(w=(0,1,-1)\) and set
\(x_1=\epsilon v\), \(x_2=\epsilon w\), \(x_3=-\epsilon(v+w)\).
The pairings of the three increments around the cycle satisfy
\[
 x_1^{\mathsf T}Mx_2=x_2^{\mathsf T}Mx_3=x_3^{\mathsf T}Mx_1
               =-3\epsilon^2.
\]
Because \(b_\kappa/(6\Delta_\kappa)\) has a positive lower
bound on \(K\), the sum of the three preferred relative
angles has circle distance at least \(c_K\epsilon^2\) from
zero for all sufficiently small positive \(\epsilon\), even
with different parameters for the three pairs.

## Controlling the deficit for arbitrary slices

The exact paired maximum for widths \(p(x),p(y)\) will be
denoted \(E_\kappa(p(x),p(y))\). Let \(D\) be the difference
between this value and the actual pairing of ordinary slices
\(A,B\). We use two consequences of a small \(D\), both uniform
for \(\kappa\in K\) and \(x,y\) small.

First, the slices are close in \(L^1\) to a common partition
into equal arcs, provided \(D\) and the mass offsets are small
enough. To see this, take a contradicting sequence with all
these offsets and \(D\) tending to zero. Extract weak-star
limits of the two slices and a limit of the kernel parameters.
The continuous kernel pairing passes to the limit by uniform
approximation with finite sums of products. Its limiting
value is the equal-arc maximum, by testing \(P(x),P(y)\) and
using weak-star compactness. The reference equality part of
JA2 applies to the soft weak-star limits as well and forces
the same ordinary equal arcs on both sides. Convergence to
indicators upgrades weak-star convergence to \(L^1\).

Second, repair \(A\) by replacing it by the best response
\(A'\) to \(B\) with its own masses, and then replace \(B\) by
the best response \(B'\) to \(A'\). Both repaired slices are
consecutive arcs in the common orientation and satisfy
\[
    \sigma(A\ne A')+\sigma(B\ne B')\le C_K\sqrt D .              \tag{RC2}
\]
Here \(\sigma\) is normalized angular measure, and inequality
of two slices means inequality of their labels. For the
first response the scores \(U_\kappa 1_{B_i}\) are uniformly
close in \(C^1\) to the equal-arc scores after a common
rotation and possibly a common reflection. The prescribed
mass assignment has prices, obtained from a supergradient
of its concave value at the interior mass vector. Every
limiting price difference must be zero, since the limiting
response maximizes the equal-arc linear pairing at equal
masses and the only such response is the reference partition.
That uniqueness follows from its strict almost-everywhere
score-winner rule. The prices are bounded after centering,
since a price difference exceeding twice the score bound
prevents a positive-mass label from winning. At the reference
endpoints the incident prices agree by continuity from both
sides. This proves that the nearby priced scores have
exactly three simple winning crossings and that the response
is an arc partition close to the reference arcs.

Uniform transversality of those crossings gives a lower
bound on the difference of the top two priced scores by a
positive constant times the distance to the three winning
endpoints. Away from fixed neighborhoods of the endpoints,
the reference winning gap is bounded below. On the set where
the label of \(A\) differs from that of \(A'\), the linear
assignment regret is at least this winning gap. The price
terms cancel because the masses agree. The integral of
distance to three points on a set of angular measure \(h\)
is bounded below by a constant times \(h^2\), since the
measure within distance \(t\) of those points is at most
\(3t/\pi\). The improvement \(F_\kappa(A',B)-F_\kappa(A,B)\)
is at most \(D\) and at least \(c_K\sigma(A\ne A')^2\).
The same estimate for the second best response proves (RC2).
The response argument also shows that both replacements
may be taken in a fixed small neighborhood of the same
reference sectors.

Define the cyclic orientation \(s(A)\in\{-1,1\}\) by the
sign of the determinant of the first angular moment vectors
of labels \(0\) and \(1\), taking the positive orientation
of \(P(0)\) as reference. For a slice of masses \(p(x)\),
define \(\varphi(A)\) by the argument of its label \(0\)
first moment minus \(s(A)\) times the center of label \(0\)
in \(P(x)\), modulo \(2\pi\). At zero determinants or zero
first moments use arbitrary conventions. In the neighborhood
above the moments are nondegenerate. Under the perturbations
in (RC2) the orientations keep their sign and
the arguments change by at most \(C_K\sqrt D\). For \(A',B'\),
these arguments are exactly their rotation angles from
\(P(x),P(y)\) or their common reflection.

The pair \(A',B'\) is also at most \(D\) below the exact
optimum. Strong angular concavity near the unique optimal
relative angle therefore bounds the distance between their
relative angle and \(s(A)u_\kappa(x,y)\) by \(C_K\sqrt D\).
Together with (RC2) this proves, whenever \(D\) is sufficiently
small,
\[
 s(A)=s(B),\qquad
 d_{\mathbb T}\bigl(\varphi(A)-\varphi(B)-s(A)u_\kappa(x,y),0\bigr)
                  \le C_K\sqrt D .                            \tag{RC3}
\]
Here \(d_{\mathbb T}\) is circle distance modulo \(2\pi\).
The thresholds in (RC3) are fixed positive constants
depending on \(K\), independent of \(\epsilon\).

## Integrating the three-slice defect

For any three small endpoint vectors \(x_1,x_2,x_3\), put
\(\Phi=u_{\kappa_{12}}(x_1,x_2)+u_{\kappa_{23}}(x_2,x_3)
+u_{\kappa_{31}}(x_3,x_1)\). For any three slices with their
respective masses,
\[
         D_{12}+D_{23}+D_{31}\ge c_K d_{\mathbb T}(\Phi,0)^2 .
                                                               \tag{RC4}
\]
If a pair deficit exceeds the fixed threshold of (RC3), use
\(d_{\mathbb T}\le\pi\) and reduce \(c_K\) if necessary. Otherwise
(RC3) gives a common orientation. The three actual relative
arguments add to zero modulo \(2\pi\). The triangle inequality
bounds \(d_{\mathbb T}(\Phi,0)\) by \(C_K\sum\sqrt{D_{ij}}\).
Cauchy--Schwarz then proves (RC4).

In particular, for the three explicit widths following (RC1),
the sum of the deficits is at least \(c'_K\epsilon^4\).
The common-angle arcs \(P(x_i)\) give an upper bound
of order \(\epsilon^4\) for this three-pair loss, by the
uniform bound on the angular second derivatives and the
order of \(u_\kappa\) in (RC1).

For a general radial mass profile sufficiently close to equal
on a compact annulus \(H\), (RC4) gives an integrated lower
bound in terms of \(\Phi\). Denote the Gaussian radial marginal
by \(\mu\), and let \(\omega_->0\) be a lower bound on the
joint radial density relative to \(\mu\times\mu\) on \(H^2\).
The loss from the exact integrated pair envelope is at least
\[
        \frac{c_K\omega_-}{3\mu(H)}
           \int_{H^3} d_{\mathbb T}\bigl(\Phi(r,s,t),0\bigr)^2
                               \,d\mu(r)d\mu(s)d\mu(t).           \tag{RC5}
\]
Indeed, integration of the left side of (RC4) against the three
independent copies of \(\mu|_H\) gives \(3\mu(H)\) times the
product radial integral of \(D\) over \(H^2\).

Now take any fixed \(0<\rho<1\) and three disjoint measurable
sets \(E_1,E_2,E_3\) of equal positive radial Gaussian
measure \(m\) in a fixed compact annulus with positive
inner radius. Prescribe \(p(x_i)\) on \(E_i\) and equal
widths elsewhere. This is a balanced mass profile because
\(x_1+x_2+x_3=0\). Each realization by a measurable
planar partition produces three ordinary slices for almost
every choice of \(r_i\in E_i\). For these slices the parameters
\(\kappa_{ij}=\rho r_ir_j/(1-\rho^2)\) lie in a fixed compact
interval \(K\subset(0,\infty)\).

Integrating the three-slice inequality with independently
drawn \(r_i\) from the normalized radial measures on \(E_i\)
gives a lower bound \(m^2 c'_K\epsilon^4\) for the sum of
pair deficits integrated over \(E_i\times E_j\) against
product radial measure. The correlated Gaussian pair of
radii has a positive continuous density relative to the
product radial measure, bounded below on the compact
annulus product. The total integrated pair envelope thus
exceeds the noise stability of every realization by at least
\(c_{\rho,K}m^2\epsilon^4\), since the other pair deficits
are nonnegative. This is a loss relative to the pair
envelope at the prescribed profile. Its comparison with
the entire radial gain and its effectiveness for arbitrary
first-contact profiles remain open.

## The cyclic width example in the first harmonic

Use increasing angular order for the labels and center the wide arc
at zero. For the widths \((1/2,1/4,1/4)\), the vector of first
coefficients \(\pi\int e^{i\theta}1_{A_j}(\theta)\,d\theta/(2\pi)\)
is
\[
  v=\left(1,\frac{-1+i}{2},\frac{-1-i}{2}\right).
\]
The other templates have vectors \((v_3,v_1,v_2)\) and
\((v_2,v_3,v_1)\). With the directed coupling \(\sum_j v_{ij}\bar
v_{i+1,j}\), all three cyclic couplings equal \(b=-1+3i/2\), so
\[
                         b^3=\frac{368+72i}{64}.
\]
Changing an angular origin multiplies the neighboring couplings
by inverse phases and preserves the product. Reversing the Fourier
convention conjugates it. The nonzero cyclic error has absolute
angle \(d=\arctan(9/46)<\pi/2\).

For arbitrary positive edge amplitudes \(a_i\), the three errors
in the first-harmonic pairing \(\sum_i a_i\cos e_i\) sum to \(d\)
modulo \(2\pi\). Take the errors in \([-\pi,\pi]\). The circle
triangle inequality gives \(\sum_i |e_i|\ge d\), so one can
choose nonnegative errors of sum \(d\), each at most the original
absolute error. This improves any tuple with \(\sum_i |e_i|>d\).
The strict concavity of cosine on \([0,d]\) gives a unique optimum
in this branch. All its errors are positive, since transferring
a small amount from a positive error to a zero error has positive
first variation. The exact loss from optimizing the three edges
independently is
\[
 \sum_i a_i\left(1-\sqrt{1-(\lambda/a_i)^2}\right),\qquad
  \sum_i\arcsin(\lambda/a_i)=d,
 \qquad 0<\lambda<(\min_i a_i)\sin d.
\]
Equal amplitudes \(a_i=a\) give \(3a[1-\cos(d/3)]\).

This first-harmonic optimum has an exact semidefinite certificate.
Gauge the coupling matrix by the optimal phases. Its directed
edges are now \(c_i+i\lambda\), where \(c_i=a_i\cos e_i\).
Imaginary row sums cancel, and the diagonal-minus-matrix slack
is \(L(c)-i\lambda J\). Here \(L(c)\) is the real weighted
triangle Laplacian and \(J_{12}=J_{23}=J_{31}=1\), \(J=-J^{\mathsf
T}\). On the orthogonal complement of the constants, \(L(c)\)
is at least \(3(\min_i a_i)\cos d\) and \(\|iJ\|=\sqrt3\).
The displayed bound on \(\lambda\) and \(\tan d=9/46<\sqrt3\)
give a positive gap on that complement. The deterministic
angular synchronization certificate therefore proves uniqueness
of the optimizer of the semidefinite relaxation.
[Deterministic synchronization certificate](THREE_CELL_ANGULAR_SYNCHRONIZATION.md#a-certificate-for-all-harmonics)

## Exact rotation optimization for the full angular kernel

The same three templates permit an all-harmonic comparison.
Prescribe cyclic shifts of \((p,(1-p)/2,(1-p)/2)\), where
\(1/3<p\le1/2\). Put \(A=2\pi p\), \(B=\pi(1-p)\),
\(\delta=(A-B)/2\) and \(c=(A+B)/2=B+\delta\). Thus
\(0<\delta\le\pi/4\), \(B\ge2\delta\) and \(B+2\delta\le\pi\).
Write \(F_{\kappa;a,b}(y)\) for the pairing of arcs of lengths
\(a,b\) at center separation \(y\) against the circular Gaussian
kernel \(K_\kappa(t)=e^{\kappa\cos t}/I_0(\kappa)\). In the
wide-centered convention a directed pair with phase difference
\(u\) has the full profile
\[
 G_\kappa(x)=2F_{\kappa;A,B}(x)+F_{\kappa;B,B}(x-\delta),
                     \qquad x=u+c.
\]
Indeed the centers of the three arcs are \(0,c,-c\). The cyclic
shift makes the two wide–narrow separations \(u+c\), and the
narrow–narrow separation \(u-2c\equiv x-\delta\).
The three variables obey \(\sum_i x_i\equiv\delta\pmod{2\pi}\),
since \(3c\equiv\delta\).

Each \(F_{\kappa;a,b}\) is even and strictly decreases with
circular distance to zero. For a given cyclic triple \(x_i\)
write \(d_{0,i}\) and \(d_{\delta,i}\) for its distances to
zero and \(\delta\). The \(y_i\in[-\delta,\delta]\) whose
distances to both centers are no larger fill the intervals
\[
 \max(-\delta,\delta-d_{\delta,i})
           \ \le y_i\le\ \min(\delta,d_{0,i}).
\]
The intervals are nonempty by the circle triangle inequality.
The same inequality gives \(\sum_i d_{0,i}\ge\delta\) and
\(\sum_i d_{\delta,i}\ge2\delta\). After capping each summand
at \(\delta\) and \(2\delta\), respectively, both bounds still
hold. The sum of the interval lower endpoints is at most
\(\delta\), and the upper sum is at least \(\delta\). One can
therefore replace the triple by such \(y_i\) with sum exactly
\(\delta\) and weakly improve every edge. Any original angle
outside \([-\delta,\delta]\) modulo \(2\pi\) gives strict
improvement on its edge.

The arc overlap formula, for \(a+b\le2\pi\), gives
\[
 F'_{\kappa;a,b}(y)
 =-C\int_{|a-b|/2}^{(a+b)/2}
               [K_\kappa(y-v)-K_\kappa(y+v)]\,dv
\]
with a common \(C>0\). This proves strict decrease on
\((0,\pi)\) and gives, writing \(K=K_\kappa\),
\[
 -G_\kappa''(x)/C =
 4K(x-\delta)+2K(x+\delta)
       -3K(x-c)-2K(x+c)-K(x+B-\delta)>0
\]
for \(|x|\le\delta\). In fact each \(K(x\pm\delta)\)
is at least both \(K(x\pm c)\) by the distance bounds
\(2\delta\le B\) and \(B+2\delta\le\pi\), with
\(K(x-\delta)>K(x-c)\). Also \(K(x+\delta)\ge
K(x+B-\delta)\). Hence with any positive edge weights and
possibly different \(\kappa_i>0\), the exact joint rotation
optimum is the unique solution of the strictly concave
problem on \([-\delta,\delta]^3\) with \(\sum_i y_i=\delta\).
Equal edge weights and equal \(\kappa\) give \(y_i=\delta/3\)
for all \(\kappa>0\). These equal-parameter triples are a
circular comparison. Distinct Gaussian radii generally give
unequal \(\kappa_i\), which the constrained problem also
allows.

Each unconstrained edge has a unique preferred \(x_\kappa\in
(0,\delta)\). Projecting into \([0,\delta]\) decreases both
circular distances, and the derivative has opposite signs at
its endpoints. If \(m_i=\min_{|x|\le\delta}(-G_{\kappa_i}''(x))>0\)
and \(w_i>0\) are the edge weights, the joint loss satisfies
\[
 \sum_i w_iG_{\kappa_i}(x_{\kappa_i})
       -\max_{\sum x_i\equiv\delta}\sum_iw_iG_{\kappa_i}(x_i)
 \ \ge\
 \frac{(\delta-\sum_i x_{\kappa_i})^2}
                       {2\sum_i(w_i m_i)^{-1}} .
\]
This follows from the concavity bound for each edge on the
clipped interval and weighted Cauchy–Schwarz.

For common \(\kappa\) the preferred angles are compatible at
exactly one \(\kappa>0\). Write \(h=\delta/3\) and \(P=2\pi/3\).
Up to a positive normalization, \(G_\kappa'(h)\) equals
\[
 2\int_0^{2h} e^{\kappa\cos t}\,dt
 -2\int_{2h}^{4h} e^{\kappa\cos t}\,dt
 -\int_{P-4h}^P e^{\kappa\cos t}\,dt
 +2\int_P^{P+2h} e^{\kappa\cos t}\,dt .
\]
The intervals are ordered since \(h\le\pi/12\). In the variable
\(\cos t\) this is a bilateral Laplace transform of an integrable
compactly supported density with two sign changes. The
[Laplace-transform zero count](THREE_CELL_TRANSVERSE_SOURCE_OSCILLATION.md#the-count-on-a-full-reference-line)
uses differentiation and Rolle's theorem and applies to this
density. The transform vanishes at zero with negative derivative.
Indeed the first-harmonic terms are proportional to
\(2ab\cos x+b^2\cos(x-\delta)\), where \(a=\sin(A/2)\)
and \(b=\sin(B/2)\) satisfy \(a=2b\cos(B/2)>b\), and their
derivative at \(h\) is \(2b\sin h(b\cos h-a)<0\).
For large positive \(\kappa\) the first positive integral
dominates, since it contains angles strictly closer to zero
than any negative term. The zero count gives a unique simple
positive root. Thus \(x_\kappa<h\) below it and \(x_\kappa>h\)
above it, by strict concavity. The all-harmonic pair preferences
are incompatible on either side of the root even though the
first-harmonic cyclic product itself is independent of
\(\kappa\).
