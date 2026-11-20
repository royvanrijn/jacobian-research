# Regular compact-complement exclusion above the critical curve

This note proves eventual actual consecutive exclusion on every compact normalized coefficient set avoiding the two simple turning points. Its cutoff depends on the fixed support and compact set. It does not itself exclude turning-point neighborhoods or prove global eventual or all-order Strong Factorial.

Fix integers a>=1,d>=2 satisfying a>(d-1)^2/4. Write

    D=a+d, b=(d-1)/2, tau=sqrt(a-b^2)>0,
    S_*=a+b-i*tau,
    Omega(S)=S^d*(D-S)/(S-a),
    psi_w(S)=S^a*(w+S^d)*exp(-S).

Let w_*=Omega(S_*), which is in the upper half-plane, and let T={w_*,conjugate(w_*)}. The actual polynomial and own-order period are

    A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k),
    I_m(w)=integral_0^infinity psi_w(S)^m dS
          =(am)!*A_m(m^d*w)/m^(D*m+1).                   (R1)

The canonical local coalescing saddle note already proves radial Airy transversality at these turning points. That local identity is an antecedent, not a new result of this note.

**Theorem AboveCriticalRegularCompactComplement.**
For every fixed support specified above and every compact
K subset of C minus T, there is an uncomputed integer N_K>=2
such that A_N,A_(N-1) have no common complex root x with
x/N^d in K whenever N>=N_K. The original integral uses each
order's own Gamma shape. There is no uniform cutoff over
supports or compact sets.

## 1. Physical roots and the no-vertical input

The physical OPEN lower wedge is -pi/d<arg S<0. It contains one simple zero zeta=(-w)^(1/d) of psi_w. Its open-wedge stationary root count is two, but for d=2 its CLOSED lower boundary has cos(pi/d)=0 and is not a decaying end. The closed-wedge flow/compactification arguments below therefore use the following local repair.

Fix a regular upper-half-plane parameter w_0 and a sufficiently small compact neighborhood U in that half-plane. Let theta_min=inf_(w in U) arg(w)>0. Choose a fixed angle-loss sigma with0<sigma<theta_min so small that both physical critical roots and zeta remain strictly inside

    C_alpha={S: -alpha<arg S<0},
    alpha=(pi-sigma)/d<pi/2,

throughout U. This is possible by their compact interior location in the original open wedge, after shrinking U. The original positive ray is unchanged. The smaller wedge has exactly those same two roots and that one zero, because it is a subset of the original wedge and contains all of them.

For S=r*exp(-i*alpha), the normalized field has B=S*Phi_S=a-S+d*S^d/(w+S^d), and

    Im B=r*sin(alpha)+d*Im(S^d*conjugate(w))/|w+S^d|^2>0.

Indeed arg(S^d)-arg(w)=-(pi-sigma)-arg(w) lies in(-2pi,-pi). On the positive angular boundary Im B<0. Thus normalized descent is strictly inward on both sides of this fixed local wedge. Now Re S>=|S|cos(alpha)>0, and infinity is uniformly decaying throughout its CLOSED angular interval. The three genuine decaying ends are0,zeta,infinity. All compactification, finite-band, cap and complete-tail arguments below refer to C_alpha on U, not to the unshrunk closed wedge when d=2. For d>=3 the same repair is harmless. Different compact parameter neighborhoods may use different fixed alpha.

The infinity cap is uniformly stable in this shrunken wedge.
Indeed V=-1/Phi_S tends to1 as |S| tends infinity, uniformly
for w in U. Choose its radius large enough that
|V-1|<cos(alpha)/2. Then
d|S|/dt=Re(V*conjugate(S)/|S|)>=cos(alpha)/2>0.
The velocity is also uniformly bounded there. This replaces
any fixed1/3 comparison that would be insufficient when
alpha is close to pi/2. The finite-end velocity expansions
V=-S/a+O(S^2) and V=-(S-zeta)+O((S-zeta)^2) give the other
two uniformly stable caps.

Differentiation gives

    Omega'(S)=-d*S^(d-1)*[(S-a-b)^2+tau^2]/(S-a)^2.      (R2)

The only nonreal branch values are w_* and its conjugate. The critical point S_* lies in the physical wedge. For example,

    d*atan(tau/(tau^2+b*(b+1)))
       <=d/(2*sqrt(b*(b+1)))<pi,

so its physical argument is strictly between -pi/d and0. The argument of w_* is

    -d*atan(tau/(tau^2+b*(b+1)))
       +atan(tau/b)+atan(tau/(b+1)).

As a function of tau its derivative is 2*d*tau^2/(a*D)>0; it starts0 and tends pi. Thus Im w_*>0. The physical root count is two, with multiplicity, by the angular-boundary exclusions and the large-w pole and large-root count. Away from w_* these two roots are simple and can be labeled holomorphically on each small parameter disc. No global labeling through a loop around w_* is assumed.

The independently audited lemma NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md supplies the following pointwise fact. At any two DISTINCT physical stationary points S_L,S_R having the same w and equal |psi_w|,

    Re(S_R-S_L)!=0.                                    (R3)

The full Corollary NoVerticalPhysicalEquimodularPair includes a=1 by its separate endpoint argument, without dividing by a-1 or K_1. Thus this note uses that full a>=1 corollary. In the present above-critical integer regime, the only a=1 case is d=2; all other geometric and contour estimates below require only a>0.

The exact stationary envelope remains

    d*w*partial_w log(lambda_R/lambda_L)=-(S_R-S_L),
    lambda_i=psi_w(S_i).                                (R4)

All logarithms in this note are local. There is no globally fixed vertical divider in the above-critical regime.

## 2. Equal-height Stokes points are regular for the complete contours

At a regular parameter, take disjoint analytic Morse discs about the two simple critical points. Truncate each descending local arm at a fixed Morse coordinate radius. On these exits the real phase has fallen by a fixed eta>0. Shrink the parameter neighborhood so that every exit loss is at least eta and

    |log|lambda_R|-log|lambda_L||<eta/3                 (R5)

whenever the center has equal critical heights.

After a descending arm leaves its Morse disc, its height is strictly below BOTH critical heights. It therefore cannot tend to either critical point: normalized descent strictly decreases height, while each potential limiting critical height is larger than the exit height. This argument applies even when the two critical factors have identical arguments. It replaces the critical-argument separation that was useful on a non-Stokes cell.

The finite-band continuation in the locally shrunken wedge C_alpha, its inward angular boundaries, stable end caps and absence of other critical points then send every arm to an actual end. The two local index-one handles join the three initial end components into the contractible compactified C_alpha disk, so the two thimbles form the same endpoint-tree basis. Their relative identities are period identities by the vanishing small and large cap integrals. For positive integer m, the differential psi_w(S)^m dS is entire.

These arguments give, uniformly on a smaller parameter neighborhood, two complete Gaussian packets

    P_i(m,w)=lambda_i(w)^m*sqrt(2*pi/m)*c_i(w)
               *(1+O(1/m)),  c_i!=0.                   (R6)

Their nonlocal arms cannot return above the fixed exit level, so those tails have uniform strict loss. The original positive cycle has coefficients0,+1,-1 in this tree basis, and at least one coefficient is nonzero. This does not assert that both saddles are active throughout the parameter plane.

At a Stokes point of unequal heights, use the finite-sublevel detour: stop the high arms at a level strictly between the heights and complete them within the lower sublevel components. The resulting high contour and the complete low thimble form a fixed basis across the wall. If the high coefficient is nonzero its Gaussian dominates; otherwise only the low Gaussian remains. This gives a zero-free neighborhood for all sufficiently large m.

At an equal-height Stokes point both factors may be active, but(R5)-(R6) provide the complete two-packet alternative directly. Thus no assertion that equal-height Stokes points do not exist is needed.

## 3. Regular actual consecutive comparison

At a regular equal-modulus parameter, use(R3) and relabel locally so that

    Re(S_R-S_L)>=c>0.

If only one packet is active, the actual period is zero-free. If both are active, their signs affect no modulus equation. Let

    E=log(lambda_R/lambda_L), B=log|c_R/c_L|.

At an actual zero of I_m, complete(R6) gives

    m*Re E(w)+B(w)=O(1/m).                              (R7)

For a common actual root x at N,n=N-1 set

    m=N-theta, w_theta=x/m^d, 0<=theta<=1.

For large N this path stays within a prescribed compact collar. The upper zero and(R4) keep Re E=O(1/N) along the path, and

    d/dtheta[m*Re E(w_theta)+B(w_theta)]
       =-Re E-Re(S_R-S_L)+O(1/N)
       <=-c+O(1/N).                                    (R8)

Both endpoint zeros would impose(R7), a contradiction. Only the values of uniform Gaussian remainders enter; no unspecified O-term is differentiated.

If the center has unequal moduli, either the active high packet dominates or only the low packet occurs. Hence a smaller neighborhood is zero-free at all large orders. This covers every regular upper-half-plane parameter, whether or not its critical phases agree. Lower-half-plane parameters follow by conjugation.

For positive real w, Omega is strictly decreasing on(a,D) by(R2). The original positive cycle has a unique nondegenerate maximum. At w=0 the maximum of S^D exp(-S) is at D. The elementary whole-ray supremum and L1 perturbation bounds supply complex zero-free neighborhoods.

For negative real w, split the original integral at its positive zero zeta. The two real subintervals have unique simple maxima, one in(0,a) and the other>D, and give exact complete two-packet expansions on a complex neighborhood. Their real-part separation is strictly positive, so the same actual comparison(R8) applies.

Consequently the regular-contour arguments prove the following precise compact-complement statement: for every fixed support in this note and every compact K subset of C minus T, there is N_K such that no common actual root has x/N^d in K at terminal orders N>=N_K. The old Gamma upper theorem captures all upper-own-order roots in |x/N^d|<=(a+d)^d. Thus every possible sequence of large-order common roots must approach one of the two turning points. This conclusion does not itself exclude such sequences.

## 4. Exact proof dependencies and remaining scope

The individual actual period, stationary map and turning-point
geometry are the explicit formulas of the retained
COALESCING_BINOMIAL_SADDLE_BRIDGE.md, Sections1--5.
The no-vertical input is precisely the full a>=1 Corollary
NoVerticalPhysicalEquimodularPair in
NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md. This source gives
only pointwise nonvanishing at equal critical magnitudes;
it is not used as a globally fixed branch order.

The inward band, exact cap, handle and complete Gaussian
MECHANISMS are written in CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md
and the finite-sublevel wall detour in
CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md.
Their cusp-only CONCLUSIONS are not applied to this regime.
Sections1--3 above explicitly verify the physical root count,
shrunk-wedge boundary signs, three decaying ends, compact
bands, absence of other physical critical points and local
real-part separation needed by those mechanisms.
The equal-height Stokes case is proved here by a fixed local
exit loss rather than by an unproved global action-phase
classification.

Only the OLD upper theorem thm:product-annulus in
GAMMA_PRODUCT_ROOT_MODULUS.tex contains the roots. Applied
with gamma=a*m+1 to the reciprocal actual polynomial it gives

    |x|<=((a+d)*m-d+1)_d<=((a+d)*m)^d.

Thus every possible sequence of sufficiently large-order
common roots must approach one of T after upper own-order
normalization. The regular theorem does NOT exclude that
remaining possibility. A separate complete local cubic
turning-neighborhood theorem is required for global
eventual noncollision. Such a local theorem is not a
premise or a conclusion of this regular module.

No finite N_K is evaluated. This note does not give a cutoff
uniform over supports, all-order SFC, an all-index assertion
about perturbed Airy roots, or a global stationary labeling.
Finite algebraic controls authenticate the displayed
identities and source scope; they do not replace the
complete contour and Gaussian arguments. All review and
execution remain internal; no formal verification or
external review is asserted.
