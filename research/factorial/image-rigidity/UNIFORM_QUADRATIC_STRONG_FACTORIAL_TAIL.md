# A uniform quadratic large-lower-exponent Strong Factorial tail

The theorem gives a uniform all-root rectangle at every finite order by exact coefficient continuation and a zero-free boundary homotopy of the actual moment polynomials.

**Theorem UniformQuadraticStrongFactorialTail.** For every integer d>=2, every integer a>=10^24*d^2, every integer N>=2 and every characteristic-zero field, the actual consecutive factorial moments of a NONZERO ordinary binomial X^a(A+B*X^d) cannot both vanish. Equivalently its two actual rational monic moment polynomials are coprime.

The analytic input is the uniform actual-phase and whole-contour bridge UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md, Sections1--6. No global Grace strip or quartic-threshold theorem is a premise of this quadratic proof. The normalizing variable is always the original coefficient ratio of the actual ordinary moments.

## 1. Real-a auxiliary moments and their exact polynomial continuation

Fix INTEGERS m>=1,d>=2, but temporarily let a=d^2/delta^2 be positive REAL for delta>0. Define

A_(m;a,d)(Z)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*Z^(m-k).

Let S have the positive Gamma density

m^(a*m+1)/Gamma(a*m+1) * S^(a*m)*exp(-m*S), S>0,

and put V=(S-a)/sqrt(a), h_delta(v)=((1+delta*v/d)^d-1)/delta. Direct expansion of the density integral gives the EXACT identity

Q_m(z,delta)=E[(z+h_delta(V))^m]
 =A_(m;a,d)(m^d*a^d*(-1+delta*z))/(m^d*a^d*delta)^m.      (R1)

This is the ACTUAL a*m+1 shape; no single common shape is used for two different orders. The polynomial is monic of degree m in z for every delta>0.

The logarithm of the centered-scaled moment-generating function is

log E exp(t*V)=-t*sqrt(a)-(a*m+1)*log(1-t/(m*sqrt(a))).

Its cumulants therefore are

kappa_1=delta/(d*m),
kappa_k=(k-1)!*[(delta/d)^(k-2)/m^(k-1)
                              +(delta/d)^k/m^k], k>=2.  (R2)

Every coefficient is a polynomial in delta, including k=2. For a finite raw moment,

E[V^j]=j!*[t^j] exp(sum_(k=1)^j kappa_k*t^k/k!),

which is a FINITE polynomial in these cumulants. Also

h_delta(v)=sum_(j=1)^d binom(d,j)*delta^(j-1)*v^j/d^j.

For fixed finite m,d only moments up to d*m occur in(R1). Thus Q_m extends to a polynomial in delta with rational coefficients in z, retaining leading coefficient1. No uniform m limit is needed to establish this extension. It is a finite exact coefficient identity for each m,d.

At delta0, kappa_1=0,kappa_2=1/m and all remaining cumulants vanish; hence

Q_m(z,0)=E[(z+G/sqrt(m))^m]
        =i^m*m^(-m/2)*He_m(-i*sqrt(m)*z),                (R3)

where G is standard real Gaussian and He_m is the monic probabilist Hermite polynomial.

## 2. Why every needed contour estimate permits REAL a

The uniform phase is

F_delta(v,z)=d^2/delta^2*[log(1+delta*v/d)-delta*v/d]
                    +log(z+h_delta(v)).

The integer-a proof's entire-integrand shorthand is not needed along the auxiliary homotopy. Every actual complex deformation used for its bulk, outside and endpoint estimates is finite and stays within a fixed bounded v rectangle, or in the physical endpoint disk |v+z/2|<.94. There

Re(1+delta*v/d)>0.

Consequently(1+delta*v/d)^(a*m) has one holomorphic logarithm branch, even for real a. The factor(z+h_delta(v))^m is still a polynomial because d,m remain integers; its zeros cause no singularity. All local phase and magnitude estimates therefore remain valid by the same finite Cauchy deformations.

The original outer pieces remain on the real axis, where1+delta*v/d>=0. At v=-d/delta its power tends to0 because a*m>0. Its endpoint is integrable, and the right endpoint decays. The uniform derivative/envelope bounds prove both long tails, with no deformation around the left branch point. This verifies the precise real-a extension needed below; no globally entire noninteger power is assumed.

In particular the WITHIN-ORDER two-packet and endpoint bounds hold for every fixed m>=10^16 and0<delta<=10^(-12), without requiring the final adjacent-order hard-range inequality. That hard-range condition is used only later when comparing two actual ordinary moments.

## 3. All initial roots lie strictly inside the rectangle

Fix eps=10^(-7) and let

R={z: |Re z|<eps, |Im z|<21/10}.

The Hermite matrix has zero diagonal and symmetric off-diagonals sqrt(j), j=1,...,m-1. Its characteristic-polynomial recurrence is He_(j+1)=x*He_j-j*He_(j-1), verifying the polynomial in(R3). It is real symmetric and irreducible tridiagonal; an eigenvector with first coordinate0 would have all coordinates0 by the row recurrence. Thus all eigenvalues are real and simple. Its maximum absolute row sum is at most2*sqrt(m-1)<2*sqrt(m). Therefore every root of Q_m(z,0) is simple, purely imaginary and has |Im z|<2. All m roots lie in R, including m1.

## 4. No boundary zero at any point of the homotopy

We now take m>=10^16 and any d>=2. All estimates are uniform for0<delta<=10^(-12).

On the horizontal boundary Im z=+/-2.1, |Re z|<=eps, the uniform outside one-packet estimate applies. It proves Q_m(z,delta)!=0. The normalized period and Q_m differ only by known nonzero Gamma/scale factors from(R1), so their zero sets agree.

For the vertical BULK boundary Re z=+/-eps, |Im z|<=1.9, let D_delta(z)=F_delta(V_+,z)-F_delta(V_-,z), with the two selected bulk saddles. At delta0 their real actions agree on the imaginary axis. The Hermite derivative is

(D_0)_z=V_+-V_-.

On the short horizontal segment from iy to+/-eps+iy its real part exceeds.59: explicitly it is the selected square root of z^2+4, whose positive real part is at least sqrt(4-y^2) for |y|<=1.9. Hence

|Re D_0(+/-eps+iy)|>.59*eps.

Both selected roots have |v|<=1+10^(-6). At a stationary point the uniform parameter derivative obeys

|partial_delta F|
 <=R^3/[6*(1-delta*R/2)^2]
           +(R^3/2)*exp(delta*R)/(1-delta*R)<1,
R=1+10^(-6).

Stationarity cancels the root derivative, giving

|D_delta-D_0|<2*delta,
|m*Re D_delta|>10^16*(.59*10^(-7)-2*10^(-12))>2>log6.    (R4)

But an actual within-order zero would require |m*Re D_delta|<log6, by the complete packet relative errors below.01 and curvature prefactor ratio below5. Thus the vertical bulk boundary has no zero. This argument uses the actual actions in each exponent and has no m*delta-smallness hypothesis.

For the upper vertical ENDPOINT boundary Re z=+/-eps,1.9<=Im z<=2.1, the Hermite matched zeta_0 is real on z=iy and has |(zeta_0)_z+i|<.03. Therefore

|Im zeta_0(+/-eps+iy)|>.97*eps.

The uniform outer complex delta disk has radius10^(-5) and |zeta|<.49 there. Cauchy's exact transfer gives

|zeta_delta-zeta_0|<=.49*delta/(10^(-5)-delta)
 <=.49/9999999<5*10^(-8).

It follows that |Im zeta_delta|>4.7*10^(-8). At any actual period zero the uniform complete-period/Robin estimate instead gives

|Im zeta_delta|<5200/m<=5.2*10^(-13).                    (R5)

These bounds contradict. The lower endpoint boundary follows by conjugation, since Q_m has real coefficients when delta>0. The top/bottom edges cover the corners, and the bulk/endpoint pieces overlap at |Im z|=1.9. Thus EVERY point of the complete rectangle boundary is zero-free for the whole closed homotopy0<=delta<=10^(-12); delta0 was covered in Section3.

## 5. The explicit uniform ALL-ROOT strip

For each fixed m,d, Q_m is monic of degree m and has coefficients continuous, indeed polynomial, in delta at0. Its boundary values never vanish by Section4. The argument principle therefore keeps its number of zeros inside R constant throughout[0,10^(-12)]. Initially that number is m, by Section3. As its degree remains exactly m, ALL roots remain inside R, counted with multiplicity. No assertion of continued simplicity is needed.

**Lemma UniformActualMomentRootRectangle.** For every d>=2, m>=10^16 and0<delta<=10^(-12), with a=d^2/delta^2 allowed real in this lemma, every root of the ACTUAL polynomial in(R1) satisfies

|Re z|<10^(-7), |Im z|<21/10.                            (R6)

This is a finite-degree topological continuation for each m, supported by uniform boundary estimates. It is not a weak-limit or sampled-zero statement. No global Grace strip or unbounded-height localization is required for(R6).

## 6. The actual two-order comparison and the quadratic threshold

Return now to INTEGER a>=10^24*d^2. Then delta=d/sqrt(a)<=10^(-12). For every N>=2 put n=N-1. The exact large-lower criterion covers a>=d*(2N-3), or equivalently N<=a/(2d)+3/2, including equality.

In the complementary hard range n>a/(2d)+1/2>=10^24*d/2+1/2, so both actual orders exceed10^16. A common coefficient root Z would give an upper normalized root z in(R6). The lower normalized parameter is EXACTLY

z_n=z+alpha*(z-delta^(-1)), alpha=(N/n)^d-1.              (R7)

Since d/n<2delta^2 and log(N/n)<1/n,

0<alpha<exp(2delta^2)-1<3delta^2.

The same bounds apply along the coupled radial path. For |z|<=3 its real shift is at most3delta+9delta^2. Thus the upper bulk region remains in the intermediate rectangle |Re z_theta|<=2*10^(-7), |Im z_theta|<=1.91, and the upper endpoint region remains in |Re z_theta|<=10^(-5),1.8999<=Im z_theta<=2.1001. Conjugation reduces to the upper half-plane where needed.

In the bulk, the exact action identity

(z-delta^(-1))*F_z=1-v/delta at a stationary point

gives dD/dtheta=-d*(V_+-V_-)/[delta*(N-theta)]. The complete uniform packets therefore give the strict adjacent-order real-action drop below-.4*sqrt(a), which excludes lower-order cancellation at an upper-order zero.

In the endpoint region, |zeta_z+i|<1/4 along(R7) gives

Im[zeta(z_n)-zeta(z)]>=alpha/(2delta)
 >=d/(2delta*N)=sqrt(a)/(2N)>=10^12/N.                    (R8)

The two actual zero bands instead bound this modulus by

5200*(1/N+1/n)<15600/N.

This contradicts(R8). Thus no hard-order common root exists. Combining with the exact soft criterion proves rational coprimality for EVERY N>=2 at the stated threshold.

Finally, the pair A_(N;a,d),A_(n;a,d) has rational coefficients. Excluding a common complex root gives a rational Bezout identity, valid over every characteristic-zero field. Monomial cases never have a zero positive factorial moment, and the zero polynomial is excluded in the theorem. This proves the numerical uniform quadratic tail.

## Exact scope

The explicit constant is10^24. The theorem covers all gaps simultaneously in the region a>=10^24*d^2. It does not prove full ordinary binomial Strong Factorial, does not classify the remaining infinite region a<10^24*d^2, and does not assert a smaller all-a>=40000 gap-two result. The real-a family is an auxiliary Gamma-moment homotopy only; the ordinary-polynomial conclusion is stated for integer exponents.


## A genuine reverse equivalence for the exceptional lower exponents

For each integer d>=2 define

Bad_d={integer a>=0: there exists N>=2 with
       gcd_Q(A_(N;a,d),A_(N-1;a,d))!=1}.

**Corollary QuadraticBadLowerExponentEquivalence.** The universal tail at the fixed constant C=10^24 is EQUIVALENT to

for every d>=2, Bad_d is contained in
{0,1,...,10^24*d^2-1}.                                  (R9)

If the tail holds then no a>=10^24*d^2 belongs to Bad_d, which proves(R9). Conversely(R9) excludes every bad actual pair at any a above the cutoff, exactly the tail. These are both directions of the equivalence, by the definition of Bad_d; neither is an attempted converse of a one-way technical certificate. The theorem independently proves both statements.

Neither statement says that Bad_d is empty. Full ordinary binomial Strong Factorial remains separate, and no bound on a bad terminal order N below the cutoff is asserted.
