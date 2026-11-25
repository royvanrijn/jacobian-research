# Uniform eventual noncollision on compact supercritical ratio intervals

Scratch mathematical candidate for whole independent review. No canonical
source or status is changed. All cutoffs below are unevaluated.

Fix real 1/4<c_-<=c_+<infinity and put I=[c_-,c_+].
The primary new claim is the following unbounded-gap statement.

**Theorem LargeGapCompactSupercriticalRatioExclusion.** There are integers
d0(I),N0(I) such that every ordinary integer support a>=1,d>=d0(I) with
c=a/(d-1)^2 in I has coprime actual polynomials A_N,A_(N-1) for every
N>=N0(I), where

    A_m(x)=sum_(j=0)^m binom(m,j)(am+1)_(dj)x^(m-j).

Equivalently there is no sequence of bad pairs with d,N both tending to
infinity and c in I. The argument below proves this sequential statement.

**Corollary CompactSupercriticalRatioUniformEventualSFC.** One unevaluated
terminal cutoff works for ALL integer supports in I, all nonzero phase
coefficients, and all characteristic-zero fields. Indeed bounded d leaves
only finitely many integer a in I; the retained407 fixed-support eventual
theorem excludes that finite set beyond its maximum cutoff. This corollary
uses407 only for that bounded-gap fallback. The new large-gap argument does
not use407, a lower root annulus, the finite CW circle law, or any constraint
m/d small. No all-order or support-global uniform cutoff is asserted here.

## 1. Exact compactification and complete actual tails

Use exactly

    s=(d-1)/2, C=a+s, L=2C/d, kappa=(d+1)/(d-1),
    p=1+2t/d, S=C+Lt=Cp, w0=C^d*kappa.

The actual original period and normalized factor are

    I_m(w)=integral_0^infinity [S^a(w+S^d)]^m exp(-mS)dS
          =(am)! A_m(m^d*w)/m^((a+d)m+1),
    I_m(w0*q)=L[C^(a+d)exp(-C)]^m J_m(q),
    J_m(q)=integral_(-d/2)^infinity Psi_(d,c)(t,q)^m dt,
    Psi_(d,c)=p^a exp(-Lt)(kappa*q+p^d).                 (A1)

The scalar is nonzero. For integer support/order the integrand is entire.
For sigma=1/d the exact cancellation is

    a Log p-Lt
      =c(1-sigma)^2[Log(1+2sigma*t)-2sigma*t]/sigma^2
          -(1-sigma)t,
    p^d=exp[Log(1+2sigma*t)/sigma].                     (A2)

This is a jointly holomorphic family on every fixed finite t/q-polydisc
for sufficiently small sigma. At sigma=0 its value is

    Psi_beta(t,q)=exp[-t^2/(2beta)](q exp(-t)+exp(t)),
    beta=1/(4c) in a compact subset B of(0,1).           (A3)

All fixed derivatives converge uniformly. No leading exponential is
replaced by its limiting value: every Gaussian or Airy exponent below is
an ACTUAL action.

Here are the whole-tail details needed to transfer finite contours.
For real t>-d/2,r>=0,H=p^d/(kappa*r+p^d) in[0,1],

    partial_t Log Psi_(d,c)(t,r)
        =2[-s-Lt+dH]/(d*p),
    |Psi_(d,c)(t,q)|<=Psi_(d,c)(t,|q|).                 (A4)

For all large d, L/d>=c_-. Choose R>=max(1,3/c_-).
The right derivative for t>=R is at most-c_-/2; the left derivative
for -d/2<t<=-R is at least1. These follow from
2(s+1)/d<=3/2,2s/d<1 and t/(1+2t/d)>=1/2 on t>=1.
Thus both complete original real tails have bounds C Psi(+-R,Q)^m/m
for bounded |q|<=Q. They include the receding left endpoint.

On fixed-height limiting rectangles the vertical end values are at most
K exp[-2c_-R^2+R+K_y]. First choose R large, then sigma small.
The actual real tails and both fixed finite vertical connectors can
therefore be put below any chosen positive selected critical modulus.
This proof uses only c_->0 and finite c_+, not c_+<1/4.
It revalidates the tail mechanism of the compact-subcritical note rather
than applying that note's theorem beyond its hypotheses.

At q=0 the limiting whole horizontal is the Gaussian
exp[t-t^2/(2beta)], with unique simple maximum t=beta. A uniform saddle
perturbation, symmetric local Gaussian calculation and the whole tails
give a zero-free actual disk |q|<r0 for all sufficiently large d,m.

The OLD Gamma product upper bound supplies every actual root with

    |q|<=((a+d)/C)^d/kappa,
    lim_(d->infinity) ((a+d)/C)^d/kappa=exp[1/(2c)]       (A5)

uniformly on I. There is one fixed upper root disk. Together with the
q0 Gaussian chart this needs no sharp lower annulus.

## 2. Model off-circle dominance and the actual circle capture

The elementary horizontal lemma for arbitrary beta>0 is used here in its
exact proved form. If q=exp(2h+i theta), t=x+iy and X=x-h, then apart
from a constant,

    V(X)=-X^2/(2beta)-hX/beta
                      +(1/2)log[cosh(2X)+cos(2y-theta)].

For h>0 reflection favors X<0; on that side every global maximum has
z=-X>h. At a stationary point v=z-h,A=sinh(2z),B=cosh(2z),
the curvature sign is the negative of

    A(B+b)-2v(1+bB), b=cos(2y-theta).

If1+bB<=0 this is positive immediately. Otherwise v<z and b<=1
bound it below by(B+1)(sinh(2z)-2z)>0.
All crossings of the stationary equation therefore have negative
derivative, giving one favored simple GLOBAL horizontal maximum.
For h<0 reflect.

For clarity, existence of the selected saddle is not inferred just from
this real maximum statement on arbitrarily prescribed horizontals.
The exact model stationary map is

    q=Q_beta(t)=exp(2t)(beta-t)/(beta+t).                (A6)

Kabluchko arXiv2203.05533v3, Theorem2.3, says that for every beta>0 and
z in the upper half-plane, zeta-beta*tan zeta=z has one simple solution
with Im zeta>Im z. Take z=-theta/2+ih for h>0 and

    t=h+i theta/2+i zeta=i beta*tan zeta.

Then Re t<0, and this is the favored stationary point.
It is simple and locally jointly analytic in beta,q. Reflection gives
h<0. On compact beta/q sets off the circle the selected saddle stays
bounded; the q0 Gaussian supplies the removable endpoint chart.
This is the sole classical model inverse theorem used for existence,
not a transfer of the model circle law or Theorem2.4 asymptotics.

Transfer each selected finite horizontal to the actual family using(A2),
the strict local curvature, strict finite outer losses and(A4).
A symmetric local Gaussian arc with complete tails gives

    J_m(q)=lambda_actual(q)^m*m^(-1/2)
               c_actual(q)(1+O(1/m)),                  (A7)

uniformly on compact q sets off the circle. The critical value and
amplitude are actual, nonzero, and uniformly controlled. Hence every
subsequential normalized-root limit for d,m->infinity lies on |q|=1.

For a common root x at N,n=N-1 the two OWN-order arguments are

    q_N=x/(N^d*w0), q_n=x/(n^d*w0)=lambda*q_N,
    lambda=(N/n)^d.                                   (A8)

Apply the actual capture to BOTH roots, using(A5). Then their moduli
tend to1, so lambda->1 and d/N->0 follow. These are consequences, not
initial assumptions on the phase compactification.

## 3. Exact model geometry on the circle

Put b=sqrt(beta(1-beta)) and

    theta_f=2acos(sqrt(beta))-2b in(0,pi),
    t_f=-ib, q_f=exp(i theta_f).                        (A9)

The stationary map derivative is

    Q_beta'(t)
       =2exp(2t)[beta(beta-1)-t^2]/(beta+t)^2.            (A10)

The only critical values in the upper and lower q half-planes are the
conjugate unit-circle turning values q_f,conjugate(q_f).
The top boundary of -pi/2<Im t<0 maps to the real axis, except the pole;
the bottom maps to the lower q half-plane. Inverse roots at bounded
nonzero targets cannot escape horizontally or to the pole.

At q=-1 the two real simple saddles are +-x_beta, where
x_beta*tanh x_beta=beta. They are the two symmetric strict global maxima
of the real horizontal. Continue this pair on the upper circle from
theta=pi down to theta_f. There is no branch value before theta_f.

Neither saddle can cross Re t=0 there. For t=-iy,0<=y<=pi/2,

    arg Q_beta(-iy)=theta(y)
        :=2atan(y/beta)-2y,
    theta'(y)=2beta/(beta^2+y^2)-2.

The maximum of theta(y) is exactly theta_f, attained at y=b. The
continuous argument lies in(-pi,theta_f], so there is no hidden2pi
winding. Thus for theta_f<theta<=pi the continued pair has nonzero
opposite real parts. Reciprocal conjugation interchanges it, so both
saddles lie on the SAME horizontal.

For h=0 the horizontal maximum lemma remains strict at every positive
stationary point. It gives at most one positive stationary maximum; its
existence forces the center to be lower. Consequently the selected pair
are exactly the two strict global maxima of the whole horizontal, with
a strict middle gap. This includes q=-1, whose midpoint can be a phase zero.

For0<=theta<theta_f there is instead a stable imaginary saddle
t=-iy,0<=y<b, with theta=theta(y) and

    Re F_tt=(y^2-b^2)/beta^2<0.                         (A11)

Its horizontal is the unique global maximum: a second positive critical
point would require an upward crossing before a downward crossing,
contradicting the same strict-curvature lemma. At q=1 it is simply t=0.
Conjugation covers the lower circle.

Thus every circle point away from the two edges has either one complete
strict simple maximum (the gap around1), or two complete strict simple
maxima (the complementary circle arcs). All margins are uniform on
compact beta sets and closed subsets avoiding the edges.

## 4. Complete actual regular-circle exclusion

At a limiting gap point the one-packet transfer(A7) is zero-free.
On a closed two-peak chart, use analytic perturbations of the two actual
saddles and split the exact integer period at their analytic midpoint.
Finite Cauchy deformations take each half to its own horizontal, with
middle connectors below BOTH peaks. The actual real tails and fixed
vertical connectors have the losses of Section1. This gives two ACTUAL
active packets

    P_i=lambda_i(q)^m*m^(-1/2)c_i(q)(1+R_i),
    |R_i|<=C/m, Re(t_R-t_L)>=c0>0.                     (A12)

The logarithms D=Log(lambda_R/lambda_L) and
b0=Log(c_R/c_L) have |Re b0|+|b0_q|<=C.
At an actual zero, Re[mD+b0]=O(1/m). Exact stationary differentiation gives

    d*q*D_q=-L*(t_R-t_L).                              (A13)

On m=N-theta,q_theta=(N/m)^d q_N, put E=mD+b0.
Its derivative is -D-L(t_R-t_L)+(d*q/m)b0_q.
The upper-zero relation and(A13) give
|Re D|<=C(1+L)/N along the path. Since d/L is bounded on I,
Re E'<=-c0 L/2 for all large members, contradicting the two endpoint
values O(1/N). No Gaussian remainder is differentiated.

A remaining common-root sequence must therefore approach one of the two
model edges. By conjugation it is enough to treat the upper edge.
This conclusion uses complete contours on the whole actual root disk,
not a zero-density or all-index model assertion.

## 5. The upper-edge horizontal selects one genuine Airy cycle

Write F_beta=Log Psi_beta near its nonzero edge critical factor.
At t_f=-ib,q_f from(A9), its exact jets are

    F_t=F_tt=0,
    F_ttt=2ib/beta^2,
    F_tttt=-2(3-2beta)/beta^2.                          (A14)

On the horizontal t=x-ib the x-dependent real part relative to x=0 is

    V(x)=-x^2/(2beta)
       +(1/2)log[(cosh(2x)+2beta-1)/(2beta)].            (A15)

It has a unique global quartic maximum0. One direct proof uses
T_beta(x)=beta*sinh(2x)/x-cosh(2x):
the numerator of its derivative, after multiplying by x^2, is

    -beta*E0(x)-2(1-beta)*x^2*sinh(2x)<0,
    E0=(1+2x^2)sinh(2x)-2x cosh(2x)>0, x>0.

Since T_beta(0)=2beta-1, the level2beta-1 has no positive solution.
The center is strictly maximal locally by(A14) and hence globally.
Every fixed outside portion |x|>=epsilon has a strict uniform loss on
compact beta. Original whole tails and large finite vertical cuts are
also below the edge height by Section1.

The analytic equations F_t=F_tt=0 for the ACTUAL family have a unique
nearby fold (t_f,actual,q_f,actual), jointly analytic in sigma,c.
Their Jacobian is nonzero by(A14) and F_tq!=0. A uniform local cubic
normal form is

    F_actual=A(q)+z^3/3-zeta(q)z, t=V(z,q),
    zeta(q_f,actual)=0.                                (A16)

Uniformity follows by the local Banach analytic normal-form construction
in the retained simple-turn proof: the model germs vary in one compact
set with nonzero cubic and mixed derivatives; take a finite cover and
reserve smaller common parameter/physical discs. This uses a local
analytic inverse only, not a global cubic-cover claim.

Set k=(beta^2/b)^(1/3)>0. Choose the branch with

    V_z=i k at the model fold,
    z=-i*x/k+[(3-2beta)*k^2/(12beta^2)]x^2+O(x^3)       (A17)

on the model horizontal. The x^2 coefficient is positive.
Thus the two fixed small horizontal endpoints map strictly into the
RIGHT two cubic decay sectors, near arguments+-pi/2 but bent right.
They deform to the two Airy rays+-pi/3, with nonzero leading amplitude.
For all small sigma this exit selection persists, after choosing the
fixed physical endpoints FIRST. The actual cubic may have a small real
part on the literal horizontal; no local actual-horizontal maximum is
assumed. Only its fixed outside loss and its two right exit labels are
used. This proves the original cycle is ONE nonzero Ai cycle, not a
Scorer cycle or an arbitrary linear combination.

The complete outer loss persists on a fixed q-disc around the actual
edge after shrinking it. It excludes every other leading contribution
from the original positive cycle. Every exponential below retains the
actual fold action.

## 6. Zero-free log lifts and genuine complex-order partial periods

The complex-order step requires a specific guard; it does not follow
from integer entire-contour identities alone.

At the model upper edge the nearest lower phase zero is

    t_zero=i(theta_f-pi)/2=-i[b+asin(sqrt(beta))].        (A18)

All other phase zeros differ by integer multiples of pi*i.
The strip between the original real axis and the fold horizontal
Im t=-b is therefore zero-free, with a COMMON positive lower margin
asin(sqrt(beta)) on the compact beta set.

On a finite enlarged rectangle containing this strip and a small collar,
Psi_beta is bounded away from0. For actual small sigma and q near the
edge the same holds by(A2); p has positive real part there as well.
Choose the analytic logarithm which agrees with the original positive-ray
branch. On the original real tails, p>0 and
Im(kappa*q+p^d)>0 because q stays in the upper half-plane with a
uniform positive imaginary part. Those logarithms match the finite
strip branch by continuation. No e^(2pi*i*m) ambiguity is introduced.

Finite local paths and the original real tails hence define exact
partial periods holomorphic in a common small complex m-sector and
in their q-charts. The real fixed outside losses persist in that sector:
on finite contours the imaginary phases are uniformly bounded, so
shrink the sector; far real phases have Gaussian/coercive real loss
and bounded factor argument. The Gamma left endpoint is integrable
because a Re m>0. Complete-tail bounds follow from(A4).

This supplies the product collars needed below, not merely a formal
complex-m replacement in an integer polynomial.

## 7. Uniform complete cubic estimates

Let zeta=zeta(q), rho=|zeta|^(1/2), K=m*rho^3.
The following estimates hold with constants independent of large d and
c in I, in one fixed sufficiently small edge disc.

Off a fixed negative-real zeta cone and for large K, the selected Ai
cycle has one nonzero Gaussian packet, hence is zero-free.
For canonical inverse arms use the polynomial properness construction:
P_chi(u)=u^3/3-chi*u, |chi|=1, and
P_chi(u(h))=P_chi(u_i)-h. Finite complex value sets have compact
preimages; velocities are O(h^(-1/2)) near a simple saddle and
O(h^(-2/3)) at infinity. At a Stokes wall a high arm is detoured in a
fixed lower sublevel, strictly below its leading height. A finite
cover of compact angular sets gives uniform losses. This is compatible
with the classical Ai single-exponential sector |arg zeta|<=pi-epsilon;
the actual class and complete outside loss were supplied in Sections5--6.
A uniformly analytic, nonzero amplitude multiplies the selected
nonzero leading value. Choose the fixed zeta radius small enough that
the outer exp(-c*m) loss is below both local selected scales.

Every large-K actual zero consequently has

    |arg(-zeta)|<epsilon.                              (A19)

In a wider negative cone put z=-sqrt(-zeta) as the common midpoint.
Its canonical phase is -4(-zeta)^(3/2)/3, strictly below both critical
heights there. The shared connector is explicit: for P(u)=u^3/3+u
the open left value-half-plane has three unbranched proper degree-one
sheets. The sheet containing the negative real u-axis approaches BOTH
critical points+-i at their boundary values+-2i/3, via
P(iv)=i(v-v^3/3),0<=|v|<=1. On Re P=-4/3, lifted vertical paths join
both left arms to u=-1. Thus the two local paths use one common low
midpoint, and their sum is the original Ai exit pair.

Lift these paths by the actual inverse V and match their original
finite exits inside the zero-free strip. Keep the real tails unchanged.
Write I here for the normalized J_m, with the nonzero support scalars in(A1) removed. Splitting that original normalized period at V(-sqrt(-zeta),q) gives I=P_L+P_R
EXACTLY, with consistent logarithms also for complex m.
Uniform Gaussian analysis on the proper arms yields

    P_i=exp(m Phi_i)*sqrt(2pi/m)c_i(1+R_i),
    |R_i|<=C/K,
    |R_(i,m)|<=C/(mK), |R_(i,q)|<=C/(K*rho^2).          (A20)

The last two are genuine Cauchy estimates on product discs
|m'-m|<c*m,|q'-q|<c*|zeta|, with support fixed.
Analytic inverse amplitudes, their derivatives, finite low connectors
and complete outer losses have common bounds by Sections5--6.
The singular common factors cancel in the amplitude ratio. Therefore

    b=b0+Log(1+R_R)-Log(1+R_L), b0=Log(c_R/c_L),
    |Re b|<=C(rho+1/K), |b_m|<=C/(mK),
    |b_q|<=C[rho^(-1)+1/(K*rho^2)].                    (A21)

For the labels compatible with the original left/right halves,
Re(S_R-S_L)>=c*L*rho, |S_R-S_L|<=C L rho.
No technical estimate is inferred just from a no-common-zero theorem.


## 8. Exact radial transversality and growing canonical fields

At the MODEL edge, exact critical differentiation gives

    q_f*zeta_q(q_f)=i*k/(2beta),                         (A22)

with k from(A17). For verification,

    F_tq=-(beta+t)^2/[2beta^2 exp(2t)],
    zeta_q=-F_tq*V_z,
    q_f*(beta+t_f)^2/exp(2t_f)=beta^2-t_f^2=beta.

The positive imaginary coefficient is bounded above and below on B.
For actual small sigma, and a small common edge neighborhood, q*zeta_q
stays in an arbitrarily narrow cone around the positive imaginary axis,
with its modulus uniformly bounded away from0. Along the actual path
q_theta=(N/(N-theta))^d q_N the image increment consequently has

    Im(zeta_n-zeta_N)>=c*(lambda-1)>=c*d/N,
    |Re(zeta_n-zeta_N)|=o(1)*Im(zeta_n-zeta_N)            (A23)

for sequences whose endpoints approach the edge. Both endpoint targets
are in the same small edge neighborhood BEFORE any canonical cone
assignment, by(A8) and the preceding regular exclusion.
The o(1) comes from sigma->0 and the endpoints tending to their limiting
model edge, where(A22) is exactly imaginary.

Suppose K_N->infinity. The upper root has the negative-cone label(A19).
Adding an almost-positive-imaginary vector cannot cancel a vector in a
narrow negative-real cone. Therefore |zeta_n|>=c|zeta_N| and
K_n>=cK_N BEFORE labeling the lower root. The lower root now has(A19).

Both endpoints lie in that cone while their difference has(A23).
It follows that their moduli are comparable and
|Delta zeta|<=C epsilon|zeta_N|. Since zeta_q is bounded and nonzero,
|Delta q|<=C epsilon|zeta_N| as well.
The image of the straight radial q segment differs from its zeta chord
by at most C|Delta q|^2. This is o(|zeta_N|), by the uniform second
derivative bound. The whole actual path thus lies in a wider negative
cone with comparable rho,K and genuine product collars.

Use the exact remainder-corrected E=mD+b, D=Phi_R-Phi_L.
At both integer endpoint zeros Re E=0. The support-fixed identities are

    E_theta=-D-DeltaS-b_m+(d*q_theta/m)*b_q,
    d*q*D_q=-DeltaS, DeltaS=S_R-S_L.                    (A24)

By(A20)--(A21), the transported amplitude terms divided by the main
L*rho are bounded by

    C*(d/L)[1/(m*rho^2)+1/(m*K*rho^3)]
       =O(rho/K+1/K^2).                               (A25)

The m-derivative contributes at most C rho^2/(L K^2).
The upper-zero estimate is |Re D|<=C(rho+1/K)/m, and its path variation
is at most C L rho/m. Since d/L is bounded, the total relative error is

    O(1/N+rho/K+1/K^2)=o(1).

The real part of E_theta is therefore strictly negative all along the
path, while the endpoint real parts both vanish. This excludes growing
canonical edge fields. No adjacent canonical step was assumed small
before BOTH endpoint captures.

## 9. Bounded canonical fields and the derivative-free comparison

It remains H_N=N^(2/3)zeta_N bounded.
If H_n were unbounded, zeta_n would be much larger than zeta_N and
would have the almost-positive-imaginary direction(A23).
Its K_n would tend to infinity, contradicting(A19).
Thus BOTH H fields are bounded.

The complete selected local cycle has the uniform two-term expansion,
after a nonzero action/scale factor is removed,

    J_m(H)=A0[Ai(H)+m^(-1/3)c1*Ai'(H)+O(m^(-2/3))].     (A26)

Here A0 is nonzero, A0,c1 and the remainder bounds are uniform over the
compact actual edge family, and their support dependence is fixed when
m changes. The proof expands the analytic inverse amplitude V_z on
z=m^(-1/3)y through its linear term on the right Ai contour.
Quadratic and zeta terms have an integrable cubic decay majorant and
relative order m^(-2/3) for bounded H. Fixed cuts and complete original
outer tails have exp(-c*m) loss. Cauchy bounds on a larger H-disc give
the same control for the required H derivatives. This is a COMPLETE
period estimate, not a local integral substituted for the original one.

Ai has only negative real zeros, all simple. The location is classical;
simplicity also follows from uniqueness for Ai''=H Ai.
Rouche on a fixed enlarged bounded field disk exhausts all actual roots
and gives finitely many separated labels. Their physical target branches
have the common-support expansion

    q_j(m)=q_f,actual+alpha_j*m^(-2/3)
                       +beta_j*m^(-1)+O(m^(-4/3)),     (A27)

uniformly. This uses the simple-root expansion in(A26) and the analytic
inverse of zeta, not a derivative of an unspecified remainder.

From(A23), including the factor n^(2/3)/N^(2/3), one has

    Im(H_n-H_N)>=c*d*N^(-1/3)+O(1/N).                   (A28)

Both actual bounded-field roots have imaginary part tending to0.
If d*N^(-1/3) failed to tend to0 this would be impossible.
Thus that scale tends to0; the full field difference tends to0 as well,
because the principal increment is almost vertical and the scaling
correction is O(1/N). Uniform separation of the finitely many real
model zeros then forces the SAME label j at both actual orders.

But their actual radial relation has

    |q_n-q_N|=|q_N|*(lambda-1)>=c*d/N,                  (A29)

whereas(A27) at the SAME label gives

    |q_j(n)-q_j(N)|
       <=C*N^(-5/3)+C*N^(-4/3)=o(d/N).                (A30)

The remainder is merely subtracted at its two endpoints.
This contradiction excludes bounded canonical fields.

## 10. Exhaustion and exact mathematical inputs

Assume there were bad pairs with N,d both tending to infinity and c in I.
Pass to c->c_* and to a limiting normalized root q_*.
The upper root disk, q0 Gaussian and model off-circle dominance give
|q_*|=1 and the same limit at the lower actual order.
Regular gap points are one-packet zero-free.
Regular two-peak arc points are excluded by Section4.
The only remaining upper/lower edges are excluded by Sections5--9,
with conjugation treating the lower edge.
This proves LargeGapCompactSupercriticalRatioExclusion.

If no pair d0,N0 existed, choose bad pairs with d>=j,N>=j.
The sequential contradiction excludes them. Finally bounded d leaves
finitely many integer a in the ratio interval, so407 supplies their
finite maximum support cutoff and proves the stated full compact-ratio
corollary. Rational Bezout then transfers to every characteristic-zero
field. For B!=0, the exact coefficient identity is

    Lambda((X^a(A+B X^d))^m)=B^m(am)! A_m(A/B);

one-term phases are nonzero monomials.

The new large-gap argument uses the raw old Gamma UPPER bound, the
self-contained all-beta horizontal lemma, and the explicitly stated
classical inverse branch theorem Kabluchko2.3 for off-circle saddle
existence. The complete finite contours, circle-gap/arc geometry,
actual single-Ai class, zero-free log lifts, full remainders and adjacent
comparison are written above.
The analytic cubic normal-form/proper-arm techniques are the same
proved local constructions in SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md
and LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md, here with ALL their uniform
family/actual-cycle hypotheses reverified in Sections5--7.
No finite CW zero law, coefficient-relative estimate, density limit,
all-index model zero perturbation, heavy finite prefix, or new lower
root annulus is used.

Kabluchko Theorem2.4 has fixed-coupling locally uniform asymptotics OFF
the circle; Remark2.6 does not supply the circle edge analysis.
DLMF9.7.5 supplies the classical Ai sector with
|arg H|<=pi-epsilon, and9.9(i) supplies its zero location.
DLMF2.4(v) allows Scorer functions for other cycles; the right-exit
proof here is indispensable. None of these classical statements by
itself proves a Gamma contour assertion.

The theorem does not address c tending to1/4, to0, or toinfinity,
nor does it claim an effective cutoff, emptiness of a finite bad set,
all-order ordinary SFC, or a uniform cutoff on ALL supports.
Everything remains scratch pending independent whole review.

