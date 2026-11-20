# An explicit absolute Airy-contour bound in the joint Airy norm

This result supplies an explicit constant in the Airy contour
estimate. It does not establish a numerical binomial Strong Factorial
cutoff, an actual complete-period remainder constant, or a fixed-parameter
inverse domain. Those estimates are separate from the contour bound proved here.

Put phi_xi(tau)=tau^3/3-xi*tau and
M(xi)=|Ai(xi)|+|Ai'(xi)|. The Ai cycle is oriented from the lower
right decay sector to the upper right decay sector; its integral is
2*pi*i*Ai(xi), as in [DLMF9.5.4](https://dlmf.nist.gov/9.5.E4).

**Lemma ExplicitJointAiryContourEnvelope.** For every xi in C there is
an Ai contour Gamma_xi, with both ends strictly in the cubic decay
sectors, such that

    integral_Gamma_xi |exp(phi_xi(tau))| |d tau|
         <2*10^7 * M(xi).                                  (A1)

The contours below also retain the local-truncation compatibility of
the retained qualitative endpoint proof: after tau=m^(1/3)*t, their
fixed-radius exits are uniformly in cubic decay subsectors when that
physical radius is sufficiently large relative to sqrt(|zeta|) and
m^(-1/3), with xi=m^(2/3)*zeta. This is not a claim that an arbitrary
specified inverse disk already meets those conditions.

## 1. Two elementary auxiliary bounds

For 0<s and t in(0,1), e^(-t)>1-t. Euler's integral therefore gives

    Gamma(s)>integral_0^1 t^(s-1)(1-t)dt=1/[s(s+1)].

In particular Gamma(2/3)>9/10 and Gamma(1/3)>9/4. If v=3^(1/6),
then 1<v<4/3. The initial-value formulas
[DLMF9.2.5--6](https://dlmf.nist.gov/9.2#ii) imply

    Bi(0)+Bi'(0)/2 <10/9+8/27=38/27<3/2.                 (A2)

For the error multiplier chi(x) in
[DLMF9.7.3](https://dlmf.nist.gov/9.7.E3), the beta integral gives

    chi(x)=pi/[2*integral_0^(pi/2) sin(theta)^x dtheta].

Consequently 0<x<=2 implies chi(x)<=2, because sin(theta)^x is at
least sin(theta)^2. We will only use x=1 and7/6. All beta/Gamma
integrals here have strictly positive real arguments.

## 2. The entire compact disk |xi|<=4

The Airy differential equation makes every Taylor coefficient of Bi
and Bi' nonnegative: the two initial coefficients are positive, the
quadratic coefficient is zero, and c_(k+3)=c_k/[(k+3)(k+2)]. Thus
their complex moduli on |xi|<=4 are bounded by their values at4.
For real x in[0,4], set B(x)=Bi(x)+Bi'(x)/2. Then

    B'(x)=Bi'(x)+(x/2)Bi(x)<=2B(x).

Equation(A2) gives B(4)<(3/2)e^8<4500. Hence on the full complex disk

    |Bi(xi)|<4500, |Bi'(xi)|<9000.

The exact Wronskian [DLMF9.2.7](https://dlmf.nist.gov/9.2.E7) yields

    1/pi <=9000*M(xi), so M(xi)>1/36000.                  (A3)

For tau=1+i*y, direct completion of the square gives

    integral_R |exp(phi_xi(1+i*y))|dy
       =sqrt(pi)*exp(1/3-Re xi+(Im xi)^2/4).

Under |xi|<=4 the exponent is at most16/3: with u=Re xi it is
at most16/3-(u+2)^2/4. Use the COMPATIBLE bent contour instead of
leaving the ends vertical:

    tau=1+x(y)+i*y,
    x(y)=0 when |y|<=6, and x(y)=(|y|-6)/sqrt(3) otherwise.

On a bend, its real phase minus the vertical phase is

    x*(1-Re xi-y^2)+x^2+x^3/3
      <=x*[5-(8/9)y^2+|y|/sqrt(3)]<=0.

The bracket is negative at |y|=6 and decreases thereafter. The arc
length factor is at most2/sqrt(3)<7/6. Since sqrt(pi)<9/5 and
e^(16/3)<225, the WHOLE bent contour has absolute integral less than

    (7/6)*(9/5)*225=945/2<500.

Together with(A3), its ratio to M is less than18,000,000. The elementary
exponential inequalities e^8<3000 and e^(16/3)<225 are checked by
positive Taylor sums and a geometric bound on the remaining tail,
without floating-point evaluations.

## 3. Explicit primary Airy error bounds when |xi|>=4

Write r=|xi| and Xi=(2/3)xi^(3/2), using the appropriate principal
branch away from the negative axis.
[DLMF9.7.17](https://dlmf.nist.gov/9.7.E17) applies to truncation
after the leading term with n=1; sigma=1/6 for Ai and sigma=0 for
Ai'. The first neglected coefficients have magnitudes5/72 and7/72.
These are exact finite remainder bounds, not an inference from a
bare asymptotic symbol.

When |arg xi|<=5*pi/6, its multiplier is less than10 for Ai and
less than7 for Ai'. In the last angular subregion, use
|cos(arg Xi)|>=1/sqrt(2), chi(7/6),chi(1)<=2, and

    sqrt(7*pi/3)*2^(7/12)+3<10,
    2*sqrt(pi)+3<7.

The other subregions have multiplier at most3. Since |Xi|>=16/3,
the relative leading-term errors are bounded respectively by

    10*(5/72)*(3/16)=25/192<1/4,
     7*(7/72)*(3/16)=49/384<1/4.                         (A4)

The explicit bent saddle contour from the retained endpoint proof has
absolute integral at most6*r^(-1/4)*exp(-Re Xi). The leading Ai term,
with(A4), gives M(xi)>=(3/4)/(2sqrt(pi))*r^(-1/4)*exp(-Re Xi).
Their ratio is less than32.

## 4. The negative-axis sector

Write xi=-r*exp(i*theta), |theta|<=pi/6, and eta=-xi. The two rotated
arguments in the exact connection formula
[DLMF9.2.12](https://dlmf.nist.gov/9.2.E12) have arguments
pi/3+theta and -pi/3+theta. Their moduli are r and their arguments
have absolute value at mostpi/2. Therefore the n=1 multiplier in
DLMF9.7.17 is at most3. The corresponding relative errors are at most

    epsilon_A=5/128, epsilon_D=7/128.                      (A5)

Set T=(2/3)eta^(3/2), U=exp(i*pi/4)*exp(-i*T), and
V=exp(-i*pi/4)*exp(i*T). The exact connection identity and its
derivative, followed by the finite error bounds, give

    2sqrt(pi)*eta^(1/4)*Ai(xi)=U+V+e_A,
    2sqrt(pi)*eta^(-1/4)*Ai'(xi)=i*(U-V)+e_D,

where |e_A|<=epsilon_A*(|U|+|V|) and
|e_D|<=epsilon_D*(|U|+|V|). Because

    |U+V|+|U-V|>=|U|+|V|,

the joint scaled norm is at least(29/32)(|U|+|V|), in particular
greater than(3/4)(|U|+|V|). Since r>=4>1, this implies

    M(xi)>=(3/4)/(2sqrt(pi))*r^(-1/4)
                       *(exp(Re Xi)+exp(-Re Xi)).         (A6)

The sum of the two exponentials is branch independent. The explicit
upper/lower negative-sector curves in the retained endpoint proof have
total absolute integral at most50*r^(-1/4) times the same exponential
sum. Thus their ratio to M is less than800/3<300. These curves are
valid on the entire closed angular sector; no sampling of Airy zeros
or assumption that Ai alone is bounded away from zero is used.

Equations(A3)--(A6) cover the complex plane. All three ratios are
strictly less than2*10^7, proving(A1). The bent and explicit rotated
contours have the same deep outgoing cubic sectors as the prior
qualitative construction, so the stated truncation compatibility is
preserved.

## 5. What this can and cannot contribute to a numerical remainder

If a TRUE divided-amplitude derivative is bounded by2300 on the local
pieces of these compatible contours inside a proved inverse domain,
then the normalized one-integration-by-parts remainder contributes

    C_local<=2300*(2*10^7)/(2*pi)<8*10^9.

This is conditional on that domain/amplitude hypothesis. To conclude
an actual C_R<=10^10 one must additionally prove that all connecting,
truncation, endpoint and original nonlocal contributions consume at
most the remaining2*10^9 in the same joint-norm normalization. Neither
that statement nor a numerical Strong Factorial threshold is proved
here. In particular the bound does not extend an auxiliary stationary
Jacobian away from the stationary points.

## Exact classical inputs read

The actual official statements inspected are DLMF9.2.1,9.2.5--7,
9.2.12;9.5.4;9.7.1--3 and9.7.17 with its post-erratum n/sigma
provisos; [5.2.1](https://dlmf.nist.gov/5.2.E1) and
[5.12.1--2](https://dlmf.nist.gov/5.12). The two explicit contour bounds6 and50
are imported from the fully written retained endpoint account,
GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md Section4. All
additional numerical inequalities here are exact rational bounds or
elementary integrals. No historical-priority, formal-verification or
external-human-review assertion is made. Updated20November2026.

## Exact reproduction

    python3 research/factorial/image-rigidity/verify_explicit_airy_contour_envelope.py \
      --repo-root . --output /tmp/explicit_airy_contour_envelope_review.json

The mandatory source ledger reuses the retained official NIST Section9.2
and Section9.7 pages. The first-order finite errors, the joint-norm
argument and the compatible contours have their stated separate scopes.
