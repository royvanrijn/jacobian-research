# An effective large-lower-exponent theorem for the ordinary gap-two family

This theorem assembles explicit actual-moment contour estimates and a piecewise complete-period Airy remainder. It concerns a numerical large-lower-exponent range, not every lower exponent or every gap.

**Theorem EffectiveGapTwoLargeLowerExponent.** Let a>=40000000000000000 be an integer. Over every characteristic-zero field, a NONZERO ordinary binomial X^a(A+B*X^2), with A,B not both zero, cannot have two consecutive vanishing factorial moments. Equivalently, for every integer N>=2 the actual degree-N and degree-(N-1) normalized moment polynomials are coprime over Q.

The field assertion concerns the actual ordinary factorial functional L(X^j)=j!, and not a fixed-parameter Gamma surrogate. The proof below uses the exact N-dependent shapes a*N+1 and a*(N-1)+1.

## 1. Actual normalization and the complementary order ranges

For m>=1 define

A_m(Z)=sum_(k=0)^m binom(m,k)*(a*m+1)_(2k)*Z^(m-k).       (E1)

If A and B are nonzero, the m-th factorial moment of X^a(A+B*X^2), after division by the nonzero factor B^m*(a*m)!, is A_m(A/B). If either A or B is zero, every nonzero-binomial moment is plainly nonzero; the zero polynomial satisfies the standard alternative. Thus actual consecutive-moment nonvanishing is equivalent to coprimality of A_N and A_(N-1).

The retained exact large-a theorem already gives pairwise coprimality through order N whenever

a>=2*(2*N-3).                                          (E2)

It remains to consider N>a/4+3/2. Put n=N-1. Both N and n then exceed10^16. For either actual order m, rescale its own moment integral by s=(original variable)/m and set

delta=2/sqrt(a), v=(s-a)/sqrt(a),
w=a^2*(-1+delta*z).

Removing the known nonzero scale factor from that actual moment gives the phase

F_delta(v,z)=4*delta^(-2)*[log(1+delta*v/2)-delta*v/2]
                 +log(z+v+delta*v^2/4).                 (E3)

More explicitly, if I_m^orig(w)=integral_0^infinity [S^a(w+S^2)]^m*exp(-m*S)dS, then

(a*m)!*A_m(Z)=m^((a+2)*m+1)*I_m^orig(Z/m^2),
I_m^orig(w)=sqrt(a)*[a^(a+2)*delta*exp(-a)]^m*I_m(z).

Every displayed scalar is nonzero. These exact factors affect neither zeros nor the following comparison. Write the normalized F-period as I_m(z); its original physical contour is inherited from the positive real moment integral and is deformed only through the entire original integrand. Within the local physical disk the logarithm branches in(E3) are single-valued and analytic.

The retained global Grace localization, actual bulk and one-packet exclusions reduce every possible common zero in the hard range to the upper or lower endpoint window. By real-coefficient conjugation it suffices to use the upper one. Both actual parameters z=z_N and z_n belong to

|Re z_m|<=10^(-5), 1.8999<=Im z_m<=2.1001,
0<delta<=10^(-8).                                      (E4)

All endpoint representations below are uniform on(E4). The lower-order parameter is the exact radially shifted actual parameter from the retained phase-shift identity, not an independent choice. Here is the explicit radial containment. The upper common zero has |Re z|<=10^(-7) and1.9<Im z<2.1. For every intermediate actual order, alpha_theta=(1-theta/N)^(-2)-1<3*delta^2. Its real shift is below3*delta+9*delta^2, while its imaginary part is multiplied by1+alpha_theta. Thus

|Re z_theta|<10^(-7)+3*delta+9*delta^2<10^(-5),
1.8999<Im z_theta<2.1*(1+3*delta^2)<2.1001.

The complete straight radial segment therefore remains in(E4).

## 2. True regular inverse and matched amplitudes

The full thin-window inverse theorem gives

F_delta(v,z)-Phi(z,delta)=phi(t)=t^3/3-zeta(z,delta)*t,
|t|<.455, |v+z/2|<.94, |zeta|<1/8.                     (E5)

The interpolating coefficients are the actual stationary Jacobians. If h=v_t, then

h=A+B*t+(t^2-zeta)*G,
|A|>.98, |A|<2, |B|<.52.                               (E6)

For each actual order these are evaluated at that order's own z_m,delta. No equality between A(z_N) and A(z_n), or between the corresponding B values, is assumed. The coefficients in(E6) are the same matched coefficients throughout the complete representation of that order.

Cauchy on the interior .45 disk gives |h|<38 for |t|<=.425. Divided-amplitude Cauchy on the .425 circle then yields

|G|<20000, |G_t|<1300000 for |t|<=.41.                   (E7)

Indeed the exact upper bounds are respectively

(.425*38)/[(.425-.41)*(.425^2-.125)]<20000,
(.425*38)/[(.425-.41)^2*(.425^2-.125)]<1300000.

The .94 physical disk contains no zero of either logarithm denominator. The actual phase, and hence exp(m*F_delta), is holomorphic there.

## 3. The finite canonical replacement and its exact orientation

At z=x+2iD, define physical endpoints u=v+z/2 by

u_+=.4,u_-=-.4 if D<=1;
u_+=X+i*s,u_-=-X+i*s if D>=1,
s=sqrt(D^2-1), X=sqrt(.16-s^2).

The certified endpoint-coordinate theorem gives

.379<|T_+|,|T_-|<.403, Re T_+,Re T_->.01,
Im T_+<-.2, Im T_->.2.                                 (E8)

Both selected labels, including their perturbations throughout(E4), are proved by complete outward-ball Rouche certificates and shared-endpoint subdiscs. The decreasing radial connectors from T_+ and T_- reach the common circle r=.41. Their cubic phase derivative with respect to radial dilation is below-.00099.

Choose the standard oriented Airy contour with the same endpoint order; its orientation may differ from the conventional Airy orientation by one fixed sign sigma in{+1,-1}. For |xi|<=4, xi=m^(2/3)*zeta, use exactly the OLD scaled contour tau=1+max(0,|y|-6)/sqrt(3)+i*y and its proved C_A=2*10^7 joint envelope. Do not substitute a different compact contour merely because it has the same homology. For |xi|>=4 use the optimized thin positive or explicit negative contours, with C_A=300.

Truncate the chosen contour at |t|=.41. Add the two right-circle arcs and radial connectors to the physical canonical endpoints. In the negative sector also join the two incoming cuts through the left part of that circle, cancelling their common incoming infinite tails. Every finite piece lies in |t|<=.41 and has a physical lift under(E5). The lift and the actual horizontal physical segment have the same endpoints, by(E8). They are both inside the simply connected physical disk .94. Cauchy's theorem therefore gives EXACT equality of their physical integrals. It is unnecessary to assert that the image of the entire actual horizontal segment lies in the cubic inverse disk.

The standard complete basis integrals, with this same sigma for both, are

integral exp(m*phi)dt=sigma*2*pi*i*m^(-1/3)*Ai(xi),
integral t*exp(m*phi)dt=-sigma*2*pi*i*m^(-2/3)*Ai-prime(xi). (E9)

The total length of the added radial and circle connectors is below10. The two radial lengths sum to less than.062, the two right arcs each have length below pi*.41/2, and the left arc has length below2*pi*.41. These coarse bounds already sum to less than10.

## 4. Uniform references, connector losses and canonical tails

Let xi=m^(2/3)*zeta and M=|Ai(xi)|+|Ai-prime(xi)|. The Hermite zeta at z=2iD is real, independently of the chosen physical endpoint. Its uniform actual perturbation from the endpoint-location proof, equation(H7), gives

|Im zeta|<E_zeta=139199/7999920000<2*10^(-5).

For the large positive sector use the STABLE reference

S_m=exp(m*Re phi(sqrt(zeta))).

One must not replace it by the larger of both critical exponentials: Ai decays in this sector. In the large negative sector use the sum of both critical exponentials. In the compact |xi|<=4 case use their sum as well. The explicit Airy estimates give the universal comparison

S_m/M<2*10^7*m^(1/6).                                   (E10)

Indeed in the compact case S_m<=450 and M>1/36000, so the ratio is below16200000. In either large sector,

M>=(3/16)*|xi|^(-1/4)*S_m,

and |xi|<m^(2/3)/8 implies(E10). This is the comparison used for every connector, including when the stable positive critical phase is negative.

The original shifted-Hermite nonlocal proof gives two separate pieces of information. First its already-Airy-normalized contribution obeys

|R_original|<60000*sqrt(m)*exp(-33*m/50000)*M.             (E11)

Second the normalized real phase at either physical cut endpoint is at least.00066 below EACH appropriate actual selected critical phase. To see the sign issue explicitly, the limiting horizontal loss is at least.0007 below -kappa, the actual horizontal perturbation is below.00003, and either actual critical action differs from its limiting action by less than.00001. Both limiting actions are at least -kappa. Consequently the endpoint lies at least.00066 below the stable positive reference when that contour is used, and below either negative-sector critical reference otherwise. This endpoint statement, not merely the integrated estimate(E11), is used on the radial and circular connectors.

The optimized large positive cuts have loss greater than1/20000 relative to their stable reference. The large negative cuts have loss greater than1/7000 relative to their own critical reference. The compact OLD contour, at m>=10^16, has its scaled shift epsilon=m^(-1/3)<1/200000; its outer |t|>=.375 portions have cubic loss at least |t|^3/5 and saddle-relative loss greater than1/100. Its two .41 cuts are in the correct right quadrants. Thus the compact C_A=2*10^7 envelope applies to exactly the contour that is cut.

Here is the whole circular comparison. On |t|=r=.41, write p=Re t. For real zeta,

g(p)=4*p^3/3-(r^2+zeta)*p

is convex on the right half-circle, because g''=8p>=0. Restoring |Im zeta|<=2*10^(-5) costs at most2*r*2*10^(-5)=.0000164 in the maximum relative to the two actual endpoint phases. The positive-sector cut loss therefore leaves

1/20000-.0000164>1/40000.

The original radial endpoint has a larger .00066 loss; hence both right-circle connectors have the same strict1/40000 loss. In the negative sector, choosing the larger critical reference only improves both endpoint inequalities. The same right-circle comparison applies.

The negative incoming cuts have p<-8/25. On that interval the real-zeta derivative satisfies

g'(p)>=4*(8/25)^2-.41^2-1/8=233/2000>0.

Moving either cut towards the negative real circle point decreases g. The two left arcs through that point therefore have loss at least1/7000-.0000164>1/10000. This explicitly closes the left connector; it is not an assumption that both incoming branches lie in the right half-plane. All connectors consequently satisfy the single bound

sup_connector |exp(m*phi)|<=S_m*exp(-m/40000).            (E12)

The finite endpoint values in the integration by parts have the stronger .00066 loss, and therefore also satisfy(E12).

For j=0 and1 the omitted complete canonical tails obey the common conservative bound

integral_tail |t|^j*|exp(m*phi)| |dt|
 <10^15*m^(-1/3)*exp(-m/40000)*M.                        (E13)

The large positive proof gives the sharper2*10^10 constant for both j. The compact OLD-contour addendum gives2*10^6. For completeness the negative weighted estimate follows directly from its proved radial estimates, without silently extending a nonconstant original amplitude to infinity. At rho=3/8 the negative tails satisfy |dt|<=2*|d|t|| and phase loss at least |t|^3/365. On each radial piece

integral_rho^infinity r*exp(-m*r^3/365)dr
 <=exp(-m*rho^3/365)*[365/(3*m*rho)
                              +365^2/(9*m^2*rho^4)].

The bracket is below10^6/m for m>=1. Including all four pieces and the reference comparison gives a constant below10^15 in(E13), with stronger exponential exp(-m/7000). The unweighted case was already proved with constant10^11 and exp(-m/14000). Thus(E13) holds for both j and all three contour regimes. Only the polynomial amplitudes1 and t are used outside the inverse disk.

## 5. Exact first-order complete-period remainder and its full budget

Let Gamma be the finite canonical replacement from Section3, including its connectors. Its endpoints are the two certified physical endpoint coordinates. By(E5)--(E7), exact integration by parts on this ENTIRE chain gives

integral_Gamma exp(m*phi)*h
 =A*integral_Gamma exp(m*phi)+B*integral_Gamma t*exp(m*phi)
  +[G*exp(m*phi)]_endpoints/m
  -integral_Gamma G_t*exp(m*phi)/m.                       (E14)

In particular there are just the two physical endpoint terms in(E14); no artificial endpoint from splitting the core and connectors is left uncancelled. The normalized complete original period is

P_m=I_m/[sigma*2*pi*i*m^(-1/3)*exp(m*Phi)].

This uses the normalized F-period from Section1, so its known nonzero original sqrt(a)/power factor has already been removed. Equations(E9) and(E14) now give the EXACT representation

P_m=A*Ai(xi)-m^(-1/3)*B*Ai-prime(xi)+R_m.                 (E15)

The common sign sigma affects neither the zero set nor any absolute estimate. There is no undocumented factor of sqrt(a), m or2*pi in(E15).

Split the G_t integral into the optimized core and the connectors. The absolute core integral is at most

C_A*m^(-1/3)*M,
C_A=2*10^7 for |xi|<=4, C_A=300 for |xi|>=4.

After division by2*pi*m^(-1/3), its contribution to R_m is at most

(1300000*C_A)/(2*pi*m)*M.

Using2*pi>6, the compact coefficient is less than(13/3)*10^12 and the large coefficient is less than65000000.

Here is an explicit bound for EVERY other term. All quantities in the following list are already divided by2*pi*m^(-1/3)*exp(m*Phi), and multiplied by M where displayed. Use(E10), the connector length<10, and2*pi>6.

- Original nonlocal part: less than60000*sqrt(m)*exp(-33*m/50000)*M by(E11).
- A/B canonical-tail replacements: since |A|+|B|<3, less than5*10^14*exp(-m/40000)*M by(E13).
- A/B connectors: |A+B*t|<3 on.41, so their contribution is below10^9*sqrt(m)*exp(-m/40000)*M.
- G_t connectors in(E14): less than5*10^13*m^(-1/2)*exp(-m/40000)*M. The raw coefficient before upward rounding is1300000*10*(2*10^7)/6.
- Both physical finite endpoints in(E14): less than2*10^11*m^(-1/2)*exp(-m/40000)*M. The raw coefficient is2*20000*(2*10^7)/6; their actual .00066 loss is weakened to1/40000 here.

The A/B connector term is needed because(E9) replaces the optimized core by the full Ai cycle; the G_t connector term is separately needed because(E14) was applied to the entire finite chain. Neither is dropped or absorbed by an unproved amplitude extension.

For m>=1 every polynomial factor in this list is at most m and every displayed exponential loss is at least1/40000. The sum of the five rounded coefficients is below10^16. Thus the full OTHER remainder E_m satisfies

|E_m|<10^16*m*exp(-m/100000)*M.                           (E16)

No contour input remains implicit in(E16): it is the sum of the five terms above. For m>=10^16, exp(x)>=x^6/720 gives

10^16*m^2*exp(-m/100000)
 <=720*10^46/m^4<1.

Consequently(E16) is smaller than M/m. Adding the local core contribution proves the following two PIECEWISE bounds for the ACTUAL complete period:

|R_m|<=10^14*M/m if |xi|<=4,
|R_m|<=10^9*M/m if |xi|>=4.                              (E17)

The compact constant is deliberately very loose. It must not be substituted globally when checking the large-xi Robin product condition.

## 6. Effective Robin continuation at every actual period zero

At a zero of P_m, equation(E15), |A|>.98 and(E17) give

|Ai(xi)/Ai-prime(xi)|<=lambda_m
 =(.52*m^(-1/3)+C/m)/(.98-C/m),                          (E18)

where C=10^14 or10^9 according to(E17). The denominator is positive. If Ai-prime were zero, the same inequality before division would force Ai=0, contradicting uniqueness for the Airy differential equation. Hence the quotient is valid at every such zero.

For |xi|>=4, use m^(2/3)>=4*10^10. The rescaled bound K=lambda_m*m^(1/3) is less than.56, so

|lambda_m^2*xi|<=.56^2/8<1/4.

For |xi|<=4, lambda_m<.011, so

|lambda_m^2*xi|<=4*.011^2<1/4.

Both ranges satisfy lambda_m^3<1/16. Uniformly the rescaled bound obeys K<2600, because the compact worst case is at most

(.52+2500)/(.98-.01)<2600.

The exact Robin-flow lemma applies with the ACTUAL complex number lambda=-Ai(xi)/Ai-prime(xi), whose modulus is bounded above by lambda_m. Its root equation is Ai(xi)+lambda*Ai-prime(xi)=0, and its continuation equation is

dxi/dlambda=-1/(1-lambda^2*xi).

The proved product and cubic guards continue it to a negative real Ai zero within distance2*|lambda|. Therefore every actual period zero in(E4) satisfies

|Im zeta|<=2*K/m<5200/m.                                (E19)

This uniform zero band includes arbitrarily large m; there is no assumption m*delta small or m bounded in terms of a.

## 7. The actual adjacent-order radial drift excludes a common zero

For the same ORIGINAL polynomial variable Z, w_N=Z/N^2 and w_n=Z/n^2. Set kappa=n/N and alpha=kappa^(-2)-1. The exact lower-order argument is

z_n=z_N+alpha*(z_N-delta^(-1)).                          (E20)

The retained endpoint reduction proves that this whole radial segment remains in(E4). Its explicit derivative bound |zeta_z+i|<1/4 gives

Im[zeta(z_n)-zeta(z_N)]>sqrt(a)/(2*N).                    (E21)

For clarity this follows directly by integrating along(E20): alpha>=2/N, |z|<=3 and delta<=1/15 imply

Im Delta zeta
 >=alpha*[(delta^(-1)-Re z)-(1/4)*|z-delta^(-1)|]
 >=alpha/(2*delta)>=sqrt(a)/(2*N).

Suppose now that the actual moments at N and n=N-1 both vanished. Both m values exceed10^16, so(E19) applies at their own actual parameters. It gives

|Im[zeta(z_n)-zeta(z_N)]|
 <2*2600*(1/N+1/n)<15600/N.

But a>=4*10^16 makes the lower bound in(E21) at least10^8/N, a contradiction. Thus the endpoint actual pair cannot have a common zero in the hard range. The bulk/outside reductions cover its complement, and the exact soft criterion(E2) covers all remaining orders. This proves coprimality of the actual pair(E1) for every N>=2.

## 8. Every characteristic-zero field and exact remaining scope

The polynomials(E1) lie in Q[Z]. The proof excludes every common complex root, so their Q-gcd is1. Their rational Bezout identity remains valid after extension to any characteristic-zero field. Rescaling back to the ordinary factorial moments proves EffectiveGapTwoLargeLowerExponent over every such field, with the nonzero-binomial condition stated in the theorem.

This theorem closes the explicit a>=4*10^16 gap-two endpoint obligation. It does not prove gap-two coprimality for every a>=40000 or all lower exponents beneath the proved threshold, does not prove the all-gap Strong Factorial Conjecture, and does not make the fixed-gap qualitative thresholds uniform in d. The proposed single complete-period constant C_R=10^10 over the whole thin window is also not established: the proof uses the stated piecewise10^14/10^9 constants. No finite scan or sampled root calculation is used to infer any of its all-order assertions.

## 9. Exact supporting inputs and reproduction scope

The following named accounts supply the separate analytic inputs used above:

- GAMMA_PRODUCT_ROOT_MODULUS.tex, cor:linear-lower-exponent: the exact soft-order criterion, including its equality boundary.
- GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md and GAP_TWO_OUTSIDE_ENDPOINT_ONE_PACKET.md: actual same-polynomial endpoint reduction and radial-segment containment at the explicit threshold.
- QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md and GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md: matched data and the true regular inverse across the full thin rectangle.
- GAP_TWO_HORIZONTAL_ENDPOINT_CUBIC_LOCATION.md: complete actual endpoint labeling, locations and decreasing radial connectors.
- EXPLICIT_AIRY_CONTOUR_ENVELOPE.md: the OLD compact C_A=2*10^7 envelope, large C_A=300 envelope, and the joint-norm lower comparisons.
- COMPACT_AIRY_PHYSICAL_CUT_ADDENDUM.md: compatible OLD compact contour cuts and weighted outer tails at the actual large orders.
- POSITIVE_AIRY_THIN_OUTER_TAIL_AND_CUT_GEOMETRY.md and NEGATIVE_AIRY_CURVE_OUTER_TAIL.md: both large-sector canonical cuts, right/left circle geometry and canonical tails. The elementary weighted-negative estimate is also displayed in Section4.
- ACTUAL_SHIFTED_HERMITE_NONLOCAL_BOUND.md: the exact entire original contour and already-normalized original nonlocal estimate, including its separate endpoint/reference loss.
- GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md: the exact Robin continuation lemma and original Airy normalization context. The present note supplies the new numerical bounds rather than importing an unspecified compactness constant.

The guarded exact rational calculations in verify_gap_two_effective_large_lower_exponent.py check the stated amplitude, connector, aggregation, Robin and drift margins, plus universal algebraic normalization identities and deliberately false alternative budgets. Those controls support the written proofs and source dependencies. They do not replace an analytic contour argument or derive all-order claims from a bounded moment prefix.
