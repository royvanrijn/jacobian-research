# Within-order simplicity in the actual central bulk

This note constructs holomorphic exact bulk periods before differentiating
their ratio. The source estimates are Sections1 and5 of
UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md. The endpoint simplicity
input used only in the final corollary is the separate scratch note
UNIFORM_ACTUAL_ENDPOINT_ROOT_SIMPLICITY.md.

**Theorem UniformActualBulkRootSimplicity.** Let a>=1,d>=2,m>=10^16 be
integers and delta=d/sqrt(a)<=10^(-10). Put

    A_m(x)=sum_(j=0)^m binom(m,j)(am+1)_(dj) x^(m-j),
    x=m^d*a^d*(-1+delta*z).

Every zero of A_m with |Re z|<=10^(-7), |Im z|<=1.9 is simple.
The statement concerns the actual degree-m polynomial and its own shape
am+1, not an adjacent-order pair.

## 1. The source domain and selected simple saddles

Use the larger OPEN bulk box

    B={|Re z|<2*10^(-7), |Im z|<1.91}.                    (B1)

For z=x+iy in B the source selected saddles V_+,V_- satisfy

    |V_i-v_0,i|<r=10^(-6),
    v_0,+/-=+/-sqrt(1-y^2/4)-iy/2,
    Re(V_+-V_-)>.59, |V_+-V_-|<3.                       (B2)

The branches are locally holomorphic in z; their labels are fixed by
their positive/negative real parts. At either saddle the curvature
C_i=-F_vv(V_i,z) has

    Re C_i>.17, |C_i|<3.                                (B3)

On local horizontal radius .05 the source supplies its Gaussian
estimate, with horizontal decay at least3u^2/50. The actual integrand
in the v coordinate is

    K_m(v,z)=exp(-m*d*v/delta)
                 *(1+delta*v/d)^(am)*(z+h(v))^m,
    h(v)=((1+delta*v/d)^d-1)/delta.

It is entire in v and polynomial in z. Integer a,m are used here.
The finite local phase F has exp(mF)=K_m wherever its chosen logarithms
are defined, and its saddle values F_i may be chosen holomorphically
on every parameter disk considered below.

## 2. Why the source's stronger-packet assignment needs a refinement

The source first cancels the common central vertical, then assigns the
remaining short connector to the stronger saddle. That assignment alone
does not define two holomorphic functions of z: the stronger saddle
depends on real parts.

Instead define the analytic midpoint

    M(z)=(V_+(z)+V_-(z))/2.

It obeys |Re M|<r and |Im M+y/2|<r. Define the EXACT partial periods

    J_-(z)=int_(-d/delta)^M(z) K_m(v,z)dv,
    J_+(z)=int_M(z)^infinity K_m(v,z)dv.                 (B4)

The tails remain on the original real line; only a finite piece
connects to M. For fixed a,d,m an entire primitive supplies the finite
piece, while the original real tails converge locally uniformly in z.
Thus both J_i are holomorphic and J_-+J_+ is exactly the original
normalized period. No logarithm is required at the midpoint, which
can itself be a zero of K_m.

For estimates deform the left path through the actual saddle height
to i*Im V_-, then to i*Im M, then horizontally to M. Reverse these two
central pieces for the right path, which then goes through V_+.
Finite Cauchy deformation is valid for the entire K_m. The two new
central connectors each have length less than2r, and all their
points satisfy

    |Re v|<r, |Im v+y/2|<r,
    |v-(-iy/2)|<2r, |v|<1.                             (B5)

The estimating paths need not depend holomorphically on z; the exact
functions (B4) already do.

## 3. Both midpoint connectors have a uniform loss

Here is the extra complete-contour estimate needed for the fixed
holomorphic split. Let v_c=-iy/2 and compare with the Hermite center
at parameter iy. The source gives its center barrier relative to
either common Hermite saddle peak S0:

    S0-Re F0(v_c,iy)>1/600,    S0=-1/2+y^2/4.            (B6)

For a point in (B5), the universal majorants in source (U2)--(U3) give

    |U(v)+v_c^2/2|<2r+delta<3*10^(-6),
    |z+h(v)-(iy+v_c)|
       <=2*10^(-7)+2r+(delta/2)*exp(delta)
        <b=1/400000.                                   (B7)

If |iy+v_c|>=.1, the change in the logarithmic real part is at most
b/(.1-b). Consequently the full center phase change is less than

    3*10^(-6)+b/(.1-b)<1/25000.                          (B8)

The same universal majorants at V_i, using |V_i-v_0,i|<r and
|iy+v_0,i|=1, give

    |Re F_i-S0|<1/50000.                                (B9)

For example the denominator perturbation is below1.3*10^(-6), its
logarithmic error is below2*10^(-6), and the quadratic plus U-tail
error is below2*10^(-6). This is also the source's stated
stationary-value comparison.

Combining (B6)--(B9), every center connector point loses more than

    1/600-1/25000-1/50000>1/1000                         (B10)

relative to BOTH actual saddle peaks.

If |iy+v_c|<.1, then |v|<.100002, Re U(v)<.01 and
|z+h(v)|<.101. Since exp(2)<8, this phase is below-1.5.
Each actual saddle peak is above-1, so the same loss (B10) follows.
A zero of the integrand has phase minus infinity and causes no problem.

This proves the missing both-packet bound. It does not infer it from
the source's stronger-packet assignment.

## 4. Exact holomorphic complete packet errors

The source local Gaussian, remaining half-horizontal, far vertical
and original real-tail estimates apply to (B4). The new center
connectors have length below2r and loss (B10). Their bound is
therefore smaller than the old displayed connector budget

    2*sqrt(m/2)*exp(-m/1000).

The complete relative errors are bounded by the sum of

    [120+240*exp(-m/5000)]/sqrt(m),
    3*sqrt(m/2)*exp(-m/10000),
    4*sqrt(m/2)*exp(-m),
    2*sqrt(m/2)*exp(-m/1000).

The source proves each is below.001 for m>10^12. Therefore throughout
B, for the PRESENT fixed holomorphic periods,

    J_i=e^(mF_i)*sqrt(2*pi/m)*C_i^(-1/2)*(1+R_i),
    |R_i|<.01.                                        (B11)

The square root is the analytic branch with positive real part,
available because Re C_i>.17. Both local horizontal Gaussians have
the same positive orientation. Thus J_i is nonzero, and R_i is
holomorphic on each parameter disk on which F_i is chosen. This is an
identity of exact periods; no estimate-only assignment defines R_i.

## 5. Within-order logarithmic derivative at an actual zero

Fix an actual zero z0 in the theorem's smaller closed box. A closed
Cauchy disk of radius

    rho_z=10^(-8)

is inside B. On its boundary |C_i|<3 and |R_i|<.01. Cauchy gives

    |(C_i)_z(z0)|<=3*10^8,
    |(R_i)_z(z0)|<10^6.                                (B12)

Set D=F_+-F_-. The exact stationary identity in source (U10), at a
single fixed order, is

    (z-delta^(-1))*D_z=-(V_+-V_-)/delta,
    D_z=(V_+-V_-)/(1-delta*z).                          (B13)

Because |z|<2 and |V_+-V_-|<3,

    |D_z-(V_+-V_-)|<6*10^(-10)/(1-2*10^(-10))
                        <10^(-8).

Hence Re D_z>.59-10^(-8)>1/2 throughout the Cauchy disk.

The logarithmic derivative difference of the two exact periods is

    E_z=(log J_+ - log J_-)_z
       =mD_z+(1/2)*[(C_-)_z/C_--(C_+)_z/C_+]
          +(R_+)_z/(1+R_+)-(R_-)_z/(1+R_-).             (B14)

No globally single-valued logarithm of the periods is required; their
logarithmic derivatives are holomorphic locally because (B11)
makes both periods nonzero. Equations (B3), (B11) and (B12) bound the
curvature contribution by

    3*10^8/(17/100)<2*10^9,

and the two remainder contributions together by

    2*10^6/(1-.01)<3*10^6.

Thus

    Re E_z>m/2-2*10^9-3*10^6>m/3>0                     (B15)

for every m>=10^16.

At a zero of J_++J_-, one has J_-=-J_+ and J_+!=0. Therefore

    (J_++J_-)_z=J_+*E_z!=0.                             (B16)

This is the required within-order simplicity, derived from actual
holomorphic periods rather than from adjacent-order coprimality.

As in the separate endpoint note, the exact nonzero normalization is

    I_m^orig(w)=Gamma(am+1)/m^((a+d)m+1)*A_m(m^d*w)
      =sqrt(a)*[a^(a+d)*delta*exp(-a)]^m*(J_-+J_+),
    w=a^d*(-1+delta*z).

Since x_z=m^d*a^d*delta!=0, (B16) proves the actual polynomial root
simple.

## 6. What the existing all-root capture actually quantifies

The canonical proof UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md,
Section2, explicitly states:

    The WITHIN-ORDER two-packet and endpoint bounds hold for every
    fixed m>=10^16 and 0<delta<=10^(-12), without requiring the final
    adjacent-order hard-range inequality.

Its Section5, Lemma UniformActualMomentRootRectangle, then states
for EVERY d>=2,m>=10^16,0<delta<=10^(-12), allowing a real in that
auxiliary lemma, that EVERY root of the actual polynomial satisfies

    |Re z|<10^(-7), |Im z|<2.1.                          (B17)

The monic finite-degree homotopy counts roots with multiplicity; it
does not assert their simplicity. Its use here is exclusively location.

In contrast, the older GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md
first restricts its adjacent comparison to hard orders
N>a/4+3/2 and uses Grace/bulk reduction for the actual common pair.
That older hard-range setup must not be inserted into the quantifiers
of the newer independent within-order capture (B17).

**Corollary UniformQuadraticHighOrderSquarefreeness.** For all integers
d>=2,a>=10^24*d^2,m>=10^16, the individual actual polynomial A_m is
squarefree over Q, and remains squarefree over every characteristic-zero
extension field.

Indeed delta<=10^(-12), so (B17) captures all roots. Those with
|Im z|<=1.9 are simple by this theorem; the upper/lower parts
1.9<=|Im z|<2.1 are simple by
UniformActualEndpointRootSimplicity. The two windows overlap.
No other root region remains. The inference is from location plus
the two independent derivative proofs, not from a homotopy assertion
of simplicity or an old soft coprimality theorem. Rational squarefreeness
gives a Bezout identity for A_m and its derivative, valid after extension
to any characteristic-zero field.

This corollary asserts high-order individual simplicity in the stated
large-lower-exponent region. It does not cover m<10^16, does not assert
all-order individual simplicity by a soft-order inference, and does
not prove support-uniform squarefreeness outside that region.

