# Consecutive exclusion in a shrinking neighborhood of a fixed simple turning point

This theorem concerns one fixed integer support and sufficiently large terminal orders. Its constants may depend on the support. It supplies no effective cutoff, estimate uniform over supports, or all-order Strong Factorial assertion.

Fix integers a>=1,d>=2 with a>(d-1)^2/4. Put D=a+d and

    psi_w(S)=S^a*(w+S^d)*exp(-S),
    I_m(w)=integral_0^infinity psi_w(S)^m dS,
    A_m(x)=sum_(k=0)^m binom(m,k)*(am+1)_(dk)*x^(m-k).

Then I_m(w)=(am)! A_m(m^d*w)/m^(Dm+1). Let S_* be the physical lower-wedge double stationary point at the upper turning value w_*. Both S_* and psi_(w_*)(S_*) are nonzero, and Phi=log psi has Phi_S=Phi_SS=0, Phi_SSS!=0, Phi_Sw!=0 there.

**Theorem FixedSimpleTurnNeighborhoodConsecutiveExclusion.**
There exist r>0 and an integer N_0, depending on this fixed support and turning point, such that for all N>=N_0 no common root x of the two actual polynomials A_N,A_(N-1) satisfies |x/N^d-w_*|<r. Conjugation would give the lower turning point as well. No effective constants or cutoff are supplied.

## 1. The local closed wedge and complete end bounds

Fix a small compact upper-half-plane neighborhood U of w_*. Use the original OPEN wedge -pi/d<arg S<0 only to identify its two stationary points with multiplicity and its one zero zeta=(-w)^(1/d). Shrink U and choose sigma with0<sigma<inf_U arg(w), so small that these selected points remain inside

    -alpha<arg S<0, alpha=(pi-sigma)/d<pi/2.

No other stationary point is introduced, since this is a subset of the original open wedge. On the lower side

    Im(S Phi_S)=r*sin(alpha)
       +d*Im(S^d*conjugate(w))/|w+S^d|^2>0,

because arg(S^d)-arg(w) is between -2pi and-pi. On the positive side it is negative. Normalized descent S_dot=-1/Phi_S is inward. Infinity is now genuinely decaying even when d=2, since cos(alpha)>0.

The three ends0,zeta,infinity have uniformly stable caps. The real phase is a log singularity at the first two and is D log|S|-Re S+O(|S|^-d) at infinity. Every closed finite phase band away from a small turning disc is compact and critical-free. The normalized field therefore gives uniform continuation through it.

For integer m, psi_w(S)^m dS is entire. Relative contour identities are exact period identities after the cap limits: the small integrals are O(r^(am+1)) and O(r^(m+1)), and the large sector arc is bounded by a polynomial times exp(-m R cos(alpha)). This argument is used only at positive integer orders when identifying the original positive-cycle period.

## 2. An exact analytic cubic normal form

There are analytic functions C(w),eta(w) and a biholomorphic coordinate S=S(t,w), on a fixed small product domain about(t,w)=(0,w_*), such that

    Phi(S(t,w),w)=C(w)+t^3/3-eta(w)*t,
    eta(w_*)=0, eta_w(w_*)!=0, G(t,w):=S_t(t,w)!=0.      (A1)

Here Phi is a single local logarithm; its domain contains no zero of psi. The two local stationary points are t=+/-sqrt(eta), and their action difference is exactly +/-4*eta^(3/2)/3.

A self-contained analytic construction is available. At w=w_* first absorb the nonzero cubic unit into t, making the critical phase exactly C_*+t^3/3. On the Banach algebra of power series with norm sum|c_j|r^j, the derivative of

    (h,C,eta) -> C+(t+h)^3/3-eta*(t+h)

at(h,eta)=(0,0) is (h,c,e)->t^2*h+c-e*t. Its inverse removes the constant and linear coefficients and divides the rest by t^2; it is bounded for each fixed r. The analytic implicit-function theorem gives h,C,eta for w near w_*. On a smaller disk Cauchy estimates make the coordinate univalent. Differentiation at the critical point and Phi_Sw!=0 give eta_w!=0. This supplies only a fixed-support local domain, as needed here.

The canonical RadialAiryParameterTransversality lemma gives, writing b=S_t(0,w_*),

    d*w_* eta_w(w_*)=b, b^3=2/Phi_SSS(S_*,w_*).          (A2)

For these actual turning points, put tau=sqrt(a-(d-1)^2/4) and S_*=a+(d-1)/2-i*tau. Then

    b^3=-i*S_*^2/tau, Im(b^3)=-Re(S_*^2)/tau<0.

The strict inequality follows by writing S_*=tau^2+c-i*tau, c=(d^2-1)/4>=3/4: Re(S_*^2)=tau^4+(2c-1)*tau^2+c^2>0. Hence b^3 is NOT real. On a ray Re eta^(3/2)=0, let q=sqrt(eta); then q^3 is purely imaginary. If Re(b*q)=0, the quotient (b*q)^3/q^3 would be real, contradicting Im(b^3)<0. Therefore the two physical critical points have

    Re(S_R-S_L)!=0 at leading order sqrt(|eta|).         (A3)

Indeed d*w*partial_w(Phi_R-Phi_L)=-(S_R-S_L) and(A1) give S_R-S_L=+/-2*d*w*eta_w*sqrt(eta). After relabeling separately on each of the three small angular cones about those rays, (A3) gives a fixed c>0 such that Re(S_R-S_L)>=c*rho, where rho=sqrt(|eta|). This is consistent with the stronger retained radial Airy angle bound; no generic Airy-combination zero-ray theorem is being imported.

## 3. The local cubic class is nonzero

At w=w_* the cubic phase has three local descending exits. Its critical-point multiplicity is two. Below the critical level there are three global end components; there is no other physical critical point anywhere in the wedge. Above that level the compactified wedge is a disk. Local resolution of the cubic point into two distinct Morse points adds exactly two index-one handles, and a critical-free exterior band retracts to the local exits and the end caps.

The three exits must reach three DISTINCT preexisting end components. If two exits reached the same component, the two handles could not turn three components into one contractible disk: one attachment would either make a loop or leave a component unjoined. Thus local cubic relative homology, with its three descending exits as endpoints, maps isomorphically to H_1(X,{0,zeta,infinity}). The original positive0-to-infinity class is nonzero. Its local cubic class is therefore nonzero as well.

This assertion is about the actual cycle, not just the dimension of an unrelated Airy model. The end assignment also identifies a path between the exits labeled0 and infinity; the resulting local contour is one of the three primitive oriented cubic Airy contours. The argument below only needs that its cubic integral is a NONZERO Airy solution.

For reviewability, the relevant elementary homology calculation is this: H_1 of a disk relative to three distinct marked endpoints is the kernel of the summation map Z^3->Z, of rank2. A local cubic disk with three descending exits has the same group. The exit-to-end bijection makes the induced map an isomorphism. The boundary of the original class is e_infinity-e_0, so it cannot be zero.

Choose the local negative exits at fixed small cubic radius r_t. At the limiting critical parameter their phases are below C_* by a fixed amount. Outside the turning disk the normalized flow has no critical point. Its three arms reach the three stable caps, with uniform strict loss along their finite connectors and tails.

The END LABELS persist uniformly, not just the count of three ends. Fix one real level strictly between the largest of the three exit heights and the critical height at w_*. Shrink the parameter disc so every chosen exit is still strictly below this level while BOTH critical heights are strictly above it. Below this level there are no critical points, so stopped inward flow identifies the sublevel with the three stable end caps throughout the parameter disc. It has exactly three components, individually labeled0,zeta(w),infinity. Each continuously moved exit remains in its original component, since its label is discrete and cannot change without crossing the chosen level or a critical point. The finite part of each connector varies continuously by compact ODE continuation; its terminal cap is uniformly stable. This proves the claimed fixed exit-to-end assignment on the entire small disc. The same construction can use the fixed physical exits from Section4, after shrinking that disc and the maximum rho once more.

Consequently the connectors persist for w in a sufficiently small complex neighborhood, and their local endpoints can be moved within fixed cubic decay sectors. Thus at INTEGER orders the original period is exactly the corresponding local cubic contour plus complete nonlocal pieces with fixed phase loss. No claim about an original noninteger-m contour identity has been made.

## 4. Individual complete thimbles and logarithmic lifts

For eta!=0 in a closed angular cone avoiding cubic Stokes directions, the two local saddles are simple. Scale

    eta=rho^2*omega, |omega|=1, t=rho*u,
    k=m*rho^3, psi_omega(u)=u^3/3-omega*u.              (A4)

The full cubic plane does NOT have compact real-phase bands, because such bands extend to infinity in oscillatory directions. The following proper-polynomial inverse construction is used instead.

Write phi_i=psi_omega(u_i), u_i=+/-sqrt(omega). A complete descending arm is a lift of

    psi_omega(u(s))=phi_i-s, s>0.                       (A4a)

At s=0 the two selected local branches resolve in sqrt(s). On a compact non-Stokes omega-set, the other critical value is never on this negative real critical-value ray. The polynomial map psi_omega is proper; hence the chosen inverse branches continue for all finite s>0 and cannot escape at a finite value of s. They cannot ramify again because the ray meets no other critical value. At infinity,

    |u(s)|~(3s)^(1/3), |du/ds|~(3s)^(-2/3),

whereas near the selected saddle |du/ds| is O(s^-1/2). These bounds are uniform on a compact non-Stokes omega-set: on each intervening finite s-interval properness bounds all roots, and separation from the other critical value bounds the inverse derivatives; the two endpoint expansions have uniform coefficients. The two chosen roots approach two distinct cubic decay ends. Thus a local Morse Gaussian and an absolutely convergent s-integral give the complete two saddle contours and uniform estimates. No compactness of an unbounded real-phase band in the cubic plane is asserted.

Truncate these curves where |t|=r_t. Since the critical phase is O(rho^3), their constant-imaginary-phase descent equations and r_t fixed give real phase at most C(w)-c_0 on those exits, for rho sufficiently small. One can see this directly from t^3/3-eta*t: on a descending end its imaginary part equals the critical imaginary part O(rho^3), and its real part is negative, hence at |t|=r_t it is at most-r_t^3/6. The exiting actual arms are below BOTH critical heights. They cannot return to either selected saddle. All remaining finite bands are uniformly regular and the global end caps are stable. This proves complete actual thimble estimates, not only estimates on the local cubic pieces.

The complete individual periods P_i can be given analytic logarithmic lifts for complex m with Re m>0. On a downward thimble choose the continuous lift satisfying

    Log psi_w(S)-Log lambda_i(w)=-s, s>=0.              (A5)

Thus the normalized exponential on that contour is exp(-m*s); its absolute decay depends on Re m, with no accumulated phase loss from winding. Local Morse Gaussian integrals are analytic on a fixed sector |arg m|<c<pi/2. At0 and zeta the local powers have integrable positive real exponents; at infinity the normalized descent coordinate has exponential decay. Truncation plus dominated convergence makes the periods holomorphic in m.

A direct inverse-ray representation makes the analytic assertions explicit. Parameterize the two arms by

    psi_w(S_(i,+/-)(s,w))=lambda_i(w)*exp(-s), s>0,

using the two inverse branches emerging in +/-sqrt(s). With an appropriate orientation,

    P_i(m,w)=lambda_i(w)^m*integral_0^infinity
           exp(-m*s)*[partial_s S_(i,+)-partial_s S_(i,-)] ds. (A5a)

There is no other physical critical value on this radial ray in the chosen non-Stokes cone. Every bounded s-interval stays in a compact physical phase band apart from the resolved starting saddle. The implicit inverse branches therefore continue holomorphically in w on a small simply connected collar. At an end their velocities are O(exp(-s/a)) at0, O(exp(-s)) at zeta, and O(1+s^C) at infinity; in fact the latter velocity tends1. The caps and their constants are uniform near the fixed turning parameter after the local turning disc is excluded. These estimates and Re m>0 justify dominated convergence and holomorphic dependence in m,w for(A5a). They also give the complete exterior Laplace tail bounds directly.

Different possible Log lambda_i lifts multiply an individual period by exp(2*pi*i*ell_i*m), where ell_i is fixed. This is absorbed into its critical logarithm. The Gaussian estimate is relative to that same lift and is unchanged. At no point is a global single-valued Log psi on the punctured wedge assumed.

The equality of the ORIGINAL positive period with an integer tree combination of P_L,P_R is needed only at positive integer orders. The entire differential and endpoint tree give coefficients0,+1,-1 there. It is unnecessary to assert that this integer homology identity is also an identity of an originally chosen principal integral at noninteger m.

## 5. Uniform singular Gaussian estimates and Cauchy collars

For a holomorphic parameter collar, rho is fixed at its positive CENTER value rho=sqrt(|eta(w_0)|), and omega(w)=eta(w)/rho^2 is allowed to vary in a small complex thickening of a closed unit omega-cone. The rescaling is therefore analytic in w; the nonholomorphic function |eta(w)| is not used as a holomorphic coordinate.

On each such closed non-Stokes omega-cone in Section4, there are constants c,C,K and rho_0>0, depending on this fixed support and cone, such that the two complete individual periods satisfy

    P_i(m,w)=exp(m*Phi_i(w))*sqrt(2*pi/m)*c_i(w)
                                  *(1+R_i(m,w)),
    |R_i|<=C/k, k=|m|*rho^3, k>=K.                     (A6)

This holds for complex m in a fixed sector and jointly on parameter discs of radius c*rho^2 about each w under consideration. All rho in a product collar are comparable.

Here are the uniformity details. In the scaled cubic coordinate the phase is EXACTLY psi_omega and the amplitude is rho*G(rho*u,w). Its ratio to its saddle value is bounded away from0 and has uniformly bounded analytic derivatives on fixed local u-discs, because G is holomorphic and nonzero on a fixed t-disc. The compact omega-parameter family has uniform local Morse maps there. A symmetric Gaussian integral eliminates its odd first correction; the analytic amplitude and fourth-order remainder give a relative O(1/k) estimate. Complete cubic tails are bounded directly from(A4a). For scaled s>=s_0 their total absolute contribution inside the chart is at most

    C*rho*|exp(k*phi_i)|*
       integral_(s_0)^infinity exp(-Re(k)*s)*(1+s)^(-2/3) ds,

where increasing C covers the bounded intervening interval. Relative to the leading rho/sqrt(|k|), this is O(exp(-c*|k|)/sqrt(|k|)), hence O(1/|k|). Extending the actual finite chart cutoff to infinity in this upper bound is legitimate and removes any apparent growing-length factor. The inverse expansion u(s)=c_j*s^(1/3)*(1+O(s^-2/3)) also makes |u(s)| strictly increasing for large s; thus each arm has a unique cut on the fixed physical circle once rho is small, in a fixed cubic decay subsector.

At the fixed physical t-exits, the phase loss is c_0>0, independent of rho. The complete exterior actual tails are uniformly bounded by a polynomial in m and rho^-1 times exp(-c_0*Re m). Relative to the leading size m^-1/2*rho^-1/2, these are uniformly O(1/k) for k>=K after K is enlarged: substitute |m|=k/rho^3; the exponential then dominates every such polynomial uniformly for0<rho<=rho_0. This is why no length or velocity factor from a growing scaled domain is simply dropped.

The holomorphic exact quotient defining R_i in(A6) is bounded on a genuine product collar

    |m'-m|<c*|m|, |w'-w|<c*rho^2.

On that collar the selected eta-sector remains non-Stokes and its scaled parameters remain in a compact cone. The preceding relative bounds apply to every point of the collar. Cauchy's formula therefore gives

    |partial_m R_i|<=C/(|m|*k),
    |partial_w R_i|<=C/(k*rho^2).                       (A7)

These are derivatives of exact holomorphic periods, not formal differentiation of a bare O-symbol.

For the leading factors, write q=sqrt(eta). Apart from a fixed unit-modulus orientation factor,

    c_R/c_L = unit*G(q,w)/G(-q,w).

The opposite Hessians contribute a fixed square-root-of-minus-one factor; their equal leading magnitudes cancel in the ratio. Hence a local logarithm b_0 of the leading ratio satisfies

    Re b_0=O(rho), partial_w b_0=O(rho^-1).             (A8)

These estimates concern the quotient; neither individual c_i is bounded at confluence.

## 6. Large scaled fields: localization and exact transport

The cubic Stokes directions have unequal real critical heights. On a closed angular region separated from Re eta^(3/2)=0, the height difference is at least c*rho^3. To justify complete uniform contours ACROSS a scaled Stokes direction, use properness in the critical-value plane, not compact cubic phase bands.

Fix one such direction and denote its lower and higher critical values by phi_L,phi_H. At the center their difference is positive real and bounded away from0. In a small closed omega-cone choose a real scaled level h_0 between their real parts, with a uniform positive margin. The negative real ray from phi_L never meets phi_H. Its two proper inverse branches therefore give the complete low contour uniformly on this cone, by(A4a).

For the high contour, lift the negative ray from phi_H only until it reaches Re value h_0. The two cutpoints are uniformly regular, since the lower critical value lies a fixed distance below this level and the high critical point has already been left by a fixed local amount. Continue in the value plane through Re value<=h_0+epsilon, where epsilon is smaller than that margin, to a negative horizontal ray tending to infinity. At the central Stokes value one may take the original ray with a fixed small semicircle around phi_L, wholly below h_0-epsilon; its lifts continue without ramification. Moving phi_L,phi_H and the initial points slightly preserves such a path, its separation from both branch values and its real-part bound. Proper polynomial inverse continuation gives uniform lifts on every finite path piece and the O(s^-2/3) velocity at infinity. A finite subcover of this small compact cone fixes all constants.

For completeness, the relevant lower-sublevel topology also follows directly from this proper cover. The preimage of the half-plane Re value<=h_0 contains exactly one simple cubic branch value. Cutting that half-plane from its branch value to its left end resolves a two-sheet component and one unramified sheet, giving two contractible components with respectively two and one cubic decay ends. The higher ramification joins these two components: the two cubic branch transpositions are distinct, since their product is the3-cycle of the roots at infinity. Hence the high cutpoints belong to different components. Selecting one end in each gives the high contour; together with the low contour it is an endpoint-tree basis. This is a proper-cover argument, rather than an assertion that a finite cubic real-phase band is compact.

All added high connectors have a fixed scaled loss below phi_H. Thus they are exponentially small compared with its local Gaussian for k large. After the lifted paths reach |t|=r_t, their negative-ray portion has the fixed physical loss established in Section4. The exterior actual finite bands, which ARE compact in the locally shrunken physical wedge, then reach the three actual ends with a fixed additional loss. The exit-to-end bijection from Section3 identifies the selected cubic ends with these actual ends. The low and high complete Gaussian contours and their endpoint-tree coefficient bounds are therefore uniform in the entire chosen Stokes cone.

The original integer period has bounded tree coefficients in that basis. If its high coefficient is nonzero, the high Gaussian dominates by exp(c*k) for k large; if not, a nonzero low coefficient gives a single nonzero Gaussian. A finite angular cover therefore shows that every actual zero with k sufficiently large lies in one of three small cones about Re eta^(3/2)=0. No global assumption of both-active saddles is used.

Use slightly wider non-Stokes cones for the analytic collars in Sections4-5. On each relabel so Re(S_R-S_L)>=c*rho as in(A3). An actual integer zero cannot have only one active Gaussian. With both coefficients nonzero, define the exact analytic correction

    b=b_0+log(1+R_R)-log(1+R_L), E=m*D+b,
    D=Phi_R-Phi_L.                                    (A9)

The logarithms exist once K is sufficiently large. The integer coefficient signs only add a constant imaginary term and do not affect the exact necessary zero equation

    Re E=0.                                           (A10)

Let x be a proposed actual common root at N,n=N-1 and put

    m=N-theta, w_theta=x/m^d, 0<=theta<=1.

This is a REAL order path. It is enough that(A10) holds at its integer endpoints. Both orders use the same individually lifted periods and the same integer tree coefficients on the small parameter collar; no noninteger positive-cycle identity is asserted.

The exact envelope d*w*D_w=-(S_R-S_L) gives

    E_theta=-D-(S_R-S_L)-b_m+(d*w/m)*b_w.               (A11)

A fixed lift difference2*pi*i*ell in D contributes only an imaginary constant along this real-order calculation. It changes neither Re E nor its derivative.

The upper zero implies Re D=O(rho/N+1/(N*k)). Along the path D_theta=-(S_R-S_L)/m keeps this bound valid with changed constants. Also(A7)-(A8) give

    b_m=O(1/(m*k)),
    b_w=O(rho^-1+1/(k*rho^2)).

Since m=k/rho^3 up to comparable constants, all real errors in(A11), other than the main negative term, are bounded by fixed multiples of

    rho^4/k + rho^3/k^2 + rho^2/k + rho/k^2.            (A12)

Choose K sufficiently large and then rho_0 sufficiently small. The right side is less than c*rho/2. Thus Re E_theta<=-c*rho/2 throughout the path, whereas both endpoint zeros would require Re E=0. This is impossible.

The collar check is essential. The relative eta-displacement in one adjacent step is O(1/(N*rho^2))=O(rho/k). First angular localization places upper zeros in narrower cones. Enlarging K and shrinking rho_0 keep the full path in the wider cone, with k and rho comparable, and within the fixed physical neighborhood. This argument excludes every common root in the large-k part of that neighborhood.

## 7. Bounded scaled fields without remainder differentiation

Put epsilon=m^-1/3 and h=m^(2/3)*eta(w). On a fixed bounded h-domain, invert eta analytically to obtain w=w(epsilon^2*h). Section3 supplies the nonzero local cubic class of the actual integer period. Its complete normalized expansion is

    J_m(h):=m^(1/3)*exp(-m*C(w))*I_m(w)
        =G_0*F(h)-epsilon*G_1*F'(h)+O(epsilon^2),        (A13)

uniformly and holomorphically in bounded complex h. Here G_0=G(0,w_*)!=0, G_1=G_t(0,w_*), and F is the nonzero canonical cubic integral of the actual local relative class.

For completeness, Taylor-expand G(epsilon*u,w)=G_0+epsilon*G_1*u+O(epsilon^2*(1+|u|^2)) on the local contour. Since w-w_*=O(epsilon^2), its parameter variation enters only at second order. The three cubic decay rays have uniformly integrable polynomial moments for bounded h. Extending from the fixed physical cut to their infinite canonical ends incurs exp(-c/epsilon^3); the complete original nonlocal pieces have the same flat bound. Thus(A13) is a full period expansion. It requires only integer m; no analyticity at epsilon=0 of an untruncated original integral is assumed.

The contour map from cubic relative homology to Airy solutions is injective. For instance, two rotated primitive contours give solutions proportional to Ai(h) and omega*Ai(omega*h), where omega is a primitive cube root of unity. Their initial-data vectors are (Ai(0),Ai'(0)) and (omega*Ai(0),omega^2*Ai'(0)); their determinant Ai(0)*Ai'(0)*(omega^2-omega) is nonzero. Hence the nonzero class in Section3 gives F not identically zero. Integration by parts in the cubic integral gives F''=h*F. Uniqueness for this ODE then shows that every zero of F is simple.

On any fixed bounded h-domain there are finitely many such zeros. Choose disjoint small discs about them and a larger boundary avoiding all model zeros. Uniform(A13) and Rouche's theorem give exactly one actual zero in each disc and no other actual zero in the prescribed smaller bounded domain. At each model zero h_j, the first correction and its uniform O(epsilon^2) remainder imply

    h_j(m)=h_j+c_j*m^-1/3+O(m^-2/3),
    w_j(m)=w_*+alpha_j*m^-2/3+beta_j*m^-1
                                      +O(m^-4/3).      (A14)

The constants and remainders are uniform over this FINITE set of roots. These are actual root locations at each integer m; no derivative of an unspecified remainder is used.

If an actual common root gives bounded upper and lower h, the two scaled fields differ by O(N^-1/3). Distinct finitely many model zeros remain separated, so the upper and lower actual roots must use the SAME label j for all sufficiently large N. But(A14) then gives

    w_j(n)-w_j(N)=O(N^-4/3)=o(1/N),

whereas the exact own-order relation requires

    w_j(n)-w_j(N)
      =[(N/n)^d-1]*w_j(N)=d*w_*/N+o(1/N).              (A15)

The turning value w_* is nonzero, so this is a contradiction.

## 8. Covering the fixed physical neighborhood

Fix the large K used in Section6. The region k<=2K has bounded h because |h|=k^(2/3), so Section7 applies with a fixed collar of that bounded model window. For sufficiently large N the lower field also belongs to that collar.

If the upper k>2K, the one-step path has comparable k>=K and lies in the wider angular and physical collars, by Section6. Angular localization and the exact ratio argument exclude it. Shrink the physical disc once to lie inside every fixed local and exterior-flow domain above. These alternatives exhaust all common roots in that fixed physical disc at sufficiently large N.

The complete thimble constructions and local exit-to-end argument in Sections3-5 establish the hypotheses used in both regimes. The alternatives above prove FixedSimpleTurnNeighborhoodConsecutiveExclusion. This is a local physical-neighborhood theorem for each fixed support. It does not itself supply a global coefficient-plane conclusion or an effective terminal cutoff.

Classical context, inspected for the specific framework: DLMF2.4(v) gives the analytic cubic change for two coalescing saddles and explicitly distinguishes contour-dependent Airy alternatives. DLMF9.5(ii), equation9.5.4, gives the canonical cubic contour; DLMF9.2 gives its ODE and nonzero initial data. The full actual-cycle, end-cap and order-coupling arguments above are additional written work, not consequences imported merely from a local CFU formula.

