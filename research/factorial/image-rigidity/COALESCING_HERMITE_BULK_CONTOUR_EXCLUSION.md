# Actual whole-contour exclusion in a uniform Hermite bulk window

This note proves an actual two-moment noncollision theorem in a specified complex coefficient window. It does not prove that every zero lies in that window. In particular the Airy endpoint and the region beyond the bulk window remain separate obligations. The large numerical constant below is chosen to make compact contour perturbations elementary.

Let d>=2 and a>=10^16*d^2 be integers, D=a+d, epsilon=1/sqrt(a), delta=d/sqrt(a)<=10^(-8), and N an integer satisfying

N>a/(2d)+3/2, n=N-1.

Thus N,n>10^12. For w in C put g(s,w)=s^a(w+s^d),

I_N(w)=integral_0^infinity g(s,w)^N exp(-Ns) ds,
J_N(w)=integral_0^infinity g(s,w)^(N-1) exp(-Ns) ds.

These are actual ordinary-factorial moment polynomials up to nonzero real factors. Write

w=a^d*(-1+delta*z), s=a+sqrt(a)*v,
h(v)=((1+epsilon*v)^d-1)/delta,
F(v,z)=epsilon^(-2)*(log(1+epsilon*v)-epsilon*v)+log(z+h(v)).

The logarithms are only local notation for absolute values and local saddle phases. The original integrand is entire in s; integer a,N ensure that contour deformations do not introduce branch cuts at its zeros or at the origin.

**Theorem ActualHermiteBulkContourExclusion.** If |Re z|<=1/20 and |Im z|<=29/20, then I_N(w) and J_N(w) do not both vanish.

The proof imports the exact local stationary branches and the exact two-order action identity in [the retained local Hermite account](COALESCING_BINOMIAL_HERMITE_BULK.md), not a presumed global saddle theorem. Sections 1--6 below establish the missing complete contour decomposition on the indicated coefficient window.

## 1. Local facts and the intermediate rectangle

The imported local lemma applies on |Re z|<=1/10, |Im z|<=3/2. Its two selected saddles V_+,V_- are within 3*delta of the roots v_+,v_- of

v^2+z*v-1=0.

At them A=-F_vv satisfies Re A>1/2 and |A|<3. On radius 1/20, |F_vvv|<6 and the horizontal real-part decay is at least u^2/5. The relative error of the local Gaussian integral is at most

[25+100*exp(-m/2000)]/sqrt(2m)

for every real m>=9999. The corresponding Gaussian leading term has modulus greater than sqrt(2/m)*exp(m*Re F(V)).

Use the intermediate rectangle

R: |Re z|<=3/40, |Im z|<=59/40.

For the actual pair, theta in [0,1], m=N-theta and kappa=m/N give z_theta by

kappa^(-d)*w=a^d*(-1+delta*z_theta).

The imported exact estimate |z_theta-z|<1/40 puts every z_theta in R. It suffices to produce complete two-packet decompositions for every z in R and every integer m>10^12.

For z=x+iy in R, Re sqrt(z^2+4)>27/20, since 4-(59/40)^2>(27/20)^2. The Hermite saddles v=t+i*eta therefore have |t|>51/80, and their moduli lie between 5/6 and 6/5. From z=1/v-v,

x=t*(1/|v|^2-1), y+eta=-eta/|v|^2.

Consequently |v|^2<=17/15. The saddle heights obey |eta|<=19*|y|/36<21/25, because |Im sqrt(z^2+4)|<=|xy|/(27/20)<=|y|/18.

## 2. A positive Hermite centre barrier

Let F_0(v,z)=-v^2/2+log(z+v), and S_0=Re F_0(v,z) at either Hermite saddle. At the centre point of the same horizontal line,

S_0-Re F_0(i*eta,z)=-t^2/2-(1/2)*log[1-t^2*(2-|v|^2)].              (B1)

If the centre is a zero of z+v, interpret the right side as positive infinity. Otherwise its logarithm has a positive argument. Put u=t^2, c=2-|v|^2. Then u>=(51/80)^2=:u_0 and c>=13/15. The function -u/2-log(1-cu)/2 increases in c. It increases in u on this region because c*(1+u)>1. The exact rational/exponential inequality

[1-(13/15)*u_0]*exp(u_0+1/50)<1

therefore proves the uniform barrier

S_0-Re F_0(i*eta,z)>1/100.                                  (B2)

The Hermite peak itself satisfies S_0>-23/25: use |v|<6/5 and log(6/5)<1/5.

On the right half-horizontal v=t+i*eta, t>=0, put Y=y+eta. The derivative of the real part of F_0 has numerator

D(t)=-t^3-2*x*t^2+(1-x^2-Y^2)*t+x,

and a nonnegative denominator. Here 1-x^2-Y^2>0. The derivative D'(t) is a quadratic with exactly one positive root. Thus this horizontal profile has at most one interior maximum, namely the selected positive-real-part saddle; any additional interior stationary point is a minimum. A possible logarithmic zero has value minus infinity and does not add a maximum. Reflection gives the same assertion for the left half-horizontal.

Near the saddle the imported local estimate at delta=epsilon=0 gives loss at least u^2/5. The centre endpoint has loss >1/100 by (B2), and the far endpoint t=3 has loss >1 (as also follows from Section 5). It follows by the preceding monotonicity that outside horizontal distance 1/25 from the saddle the Hermite phase is at most

S_0-1/3125.                                                (B3)

## 3. Uniform compact perturbation

On |v|<=16/5 the elementary Taylor estimates give

|h(v)-v|<=delta*(16/5)^2/[2*(1-(16/5)*delta)]<6*delta,
|epsilon^(-2)*(log(1+epsilon*v)-epsilon*v)+v^2/2|<6*delta.    (B4)

Compare an actual horizontal line with its Hermite line by changing the height by at most 3*delta. The quadratic part changes by less than 10*delta, and z+h(v) changes by less than 9*delta. Whenever the Hermite denominator has modulus at least 1/10, the logarithmic real part changes by less than 91*delta. Thus the total phase change is less than 120*delta. The same bound compares each actual stationary value S with its Hermite value S_0.

If the Hermite denominator has modulus below 1/10, the actual denominator has modulus below 101/1000. Since |eta|<21/25 and the height perturbation is at most 3*delta, the real part of the first term in F is below 3/8. The exact inequality

(101/1000)*exp(15/8)<1

gives Re F<-3/2 on this exceptional near-zero set. Meanwhile every actual saddle peak S>-1.

Outside actual horizontal distance 1/20 from its saddle, the corresponding Hermite distance exceeds 1/25. Equations (B3)--(B4), and

1/3125-240*10^(-8)>1/4000,

therefore show that every such actual horizontal point has phase at most S-1/4000. Each nonlocal half-horizontal has length below 3.

## 4. The short connector belongs to the stronger packet

For the Hermite heights and y>=1/10 one has

-19*y/36<=eta<=-17*y/36,
17*y/36<=y+eta<=19*y/36.

On the interval between the two heights, coarsen these to -0.53*y<=eta and 0.47*y<=y+eta<=0.53*y. The centre-line derivative satisfies

d/deta Re F_0(i*eta,z)=eta+(y+eta)/(x^2+(y+eta)^2)>0,

because the exact rational inequality

(47/100)/[(3/40)^2+((53/100)*(59/40))^2]>53/100

holds. Hence the largest centre-line phase between the saddle heights occurs at an endpoint. Equation (B2) puts it at least 1/100 below the larger Hermite peak. For negative y, conjugate the argument.

When |y|<1/10, every intermediate height has |eta|<19/360, the denominator |z+i*eta| is below 1/4, and the Gaussian term is below 1/100. Therefore the centre phase is less than -6/5, while the larger saddle peak exceeds -23/25.

The actual connector interval differs from the Hermite interval by at most 3*delta at each endpoint. Project any actual height onto the Hermite interval and apply the compact comparison in Section 3. In its small-denominator alternative the phase is below -3/2. Otherwise the loss from the larger actual peak is at least

1/100-240*delta>1/200.

The connector length is below 2. It is therefore an exponentially small error for the STRONGER packet. Nothing here claims that it is smaller than both packets.

## 5. Outer real tails and the far verticals

These estimates hold already for delta<=1/100 and epsilon<=delta/2. They are included explicitly because local saddle estimates do not control the full original contour.

For real v>=3, h(v)>=v, and Bernoulli's inequality gives

h'(v)/h(v)<= (1/v+delta)/(1+epsilon*v).

Since |z|<8/5,

Re F_v(v,z)<= [-v+(15/7)*(1/v+delta)]/(1+epsilon*v)<-2.       (B5)

For the last inequality multiply by the positive denominator; its worst bound is already negative at v=3, epsilon=1/200, delta=1/100 and decreases with v. Also

Re F(3,z)<=-900/203+log 5<-5/2.

Thus the right real tail has absolute integral at most exp(m*Re F(3,z))/(2m).

For the left real tail write v=-r, 3<=r<1/epsilon. Here h is negative, h'<=1, and

|h(-r)|>=(1-exp(-3*delta))/delta>=3-(9/2)*delta>=591/200.

It follows that Re F_v(-r,z)>=r-200/271>2. The left endpoint value satisfies Re F(-3,z)<=-9/2+log 5<-5/2, so the same integral bound holds on the left tail. At the original endpoint v=-1/epsilon the entire integrand vanishes; the logarithmic notation tends to minus infinity.

On either vertical line Re v=+3 or -3 between the real axis and an actual saddle height, |Im v|<21/25 and |v|<16/5. The elementary compact bound

Re F(v,z)<=-2592/625+4096/73800+log 5<-12/5                (B6)

uses the logarithm Taylor remainder at epsilon<=1/200 and log 5<5/3. Each vertical has length below 1. Since every peak exceeds -1, the real tails and far verticals assigned to either packet together have relative size below

4*sqrt(m/2)*exp(-m)

for m>=10^12.

## 6. Exact contour and complete error budget

Deform the entire original contour v in [-1/epsilon,infinity) to the path

real [-1/epsilon,-3]; vertical at -3 to Im V_-;
horizontal to i*Im V_-; short vertical to i*Im V_+;
horizontal to 3+i*Im V_+; vertical at 3 back to the real axis;
real [3,infinity).

This deformation cancels the common central vertical before any absolute value is taken. It is justified on finite rectangles by Cauchy's theorem; the right endpoint remains on the original real tail. Assign the short connector to whichever selected saddle has the larger real phase. Each resulting packet contains its full horizontal local Gaussian segment with positive orientation.

The local Gaussian error is at most

[25+100*exp(-m/2000)]/sqrt(2m)<1/1000.

The nonlocal half-horizontal contributes at most

3*sqrt(m/2)*exp(-m/4000)<1/1000.

The assigned real tail and far vertical contribute at most

4*sqrt(m/2)*exp(-m)<1/1000.

For the stronger packet the short connector contributes at most

2*sqrt(m/2)*exp(-m/200)<1/1000.

All four displayed inequalities hold at m=10^12 and thereafter by direct differentiation; elementary exponential-series bounds certify that endpoint. Hence each complete packet has relative error below 1/100, in particular below 1/5. Its leading Gaussian has the exact actual stationary phase and a curvature prefactor whose ratio between the packets is less than sqrt(6).

The warning about common verticals is essential. For example at z=3i/2 the Hermite phase on the original centre v=0 is log(3/2), whereas each selected saddle phase is 1/16. Thus bounding the two original real half-contours separately would include a central vertical exponentially larger than the saddle packets. The deformation above removes precisely that invalid step.

## 7. Exact two-order separation and the actual pair

The exact change of variables s=(n/N)*r gives

J_N(w)=(n/N)^(D*n+1)*I_n((N/n)^d*w).                       (B7)

For theta in [0,1], let S_theta,+ and S_theta,- be the corresponding original-s stationary branches and put

Delta_theta=Re[(N-theta)*log g(S_theta,+,w)-N*S_theta,+
              -(N-theta)*log g(S_theta,-,w)+N*S_theta,-].

The imported exact action identity says

d/dtheta [Delta_theta/(N-theta)]
    =-N*Re(S_theta,+-S_theta,-)/(N-theta)^2.

The saddle separation exceeds (6/5)*(1-theta/N)*sqrt(a). Consequently

Delta_1 <(n/N)*Delta_0-sqrt(a).                             (B8)

If I_N(w)=0, its two complete packets with errors below 1/5 and curvature ratio below sqrt(6) force |Delta_0|<log 4. Equations (B7)--(B8) make the ratio of the two lower-order packet magnitudes less than 16*exp(-sqrt(a))<1, including their relative errors. Therefore J_N(w) is nonzero. This proves ActualHermiteBulkContourExclusion.

## 8. Exact scope and the remaining global problem

The theorem concerns actual ordinary-factorial consecutive moments, not two polynomials with a common fixed Gamma shape. It supplies a whole-contour estimate, not just a local cubic or Gaussian approximation. It is uniform in the specified unbounded integer a,d,N range.

The hypothesis |Re z|<=1/20, |Im z|<=29/20 is essential to the proved account. No global zero localization is imported. The separate Grace--Szego gap-two localization gives |Re z|<1/20 and |Im z|<4 for sufficiently large a, and therefore leaves the parts 29/20<|Im z|<4 outside this theorem. The turning-point endpoints, their actual cycle selection, and quantitative adjacent-order Airy control are not solved here. No all-gap-two or full binomial SFC conclusion follows yet.


## Exact reproduction

    python3 research/factorial/image-rigidity/verify_coalescing_hermite_bulk_contour_exclusion.py \
      --repo-root . --output /tmp/coalescing_bulk_contour_review.json

The exact supporting comparisons are guarded by the proof/checker ledger
and the retained local Hermite proof/checker. They support the written
all-parameter contour proof; finite inequalities alone do not replace it.

Updated 20 November 2026.
