# Consecutive exclusion throughout a fixed neighborhood of an integer cusp

This note proves the result for each FIXED integer cusp. Constants may depend on that cusp. The retained fixed-cusp bounded-window theorem, or its stronger uniform version, supplies the central scale. The present argument treats unbounded scaled fields by complete horizontal-contour packets and derivatives of EXACT analytic amplitude functions. No differentiation of an unspecified O-term or same-index assertion is used.

Fix s>=1, a=s^2, d=2s+1, S_c=s(s+1), L=sqrt(S_c), and w_c=S_c^d*(s+1)/s. Let

    Phi(S,w)=a Log(S)+Log(w+S^d)-S,
    I_m(w)=integral_0^infinity exp(m Phi(S,w)) dS.

For w in a fixed small disc about w_c, choose the logarithms positive on the original positive ray at real w. Fix 0<epsilon_d<pi/(2d); on the narrow S-sector |arg S|<epsilon_d, Re(S^d)>0. If Re(w)>0 then w+S^d is nonzero there. Thus the branches exist throughout that sector, not throughout an unjustifiably large right half-plane. For complex m in a sufficiently small fixed sector about the positive real axis the original integral converges at0 and infinity and is holomorphic in (m,w), locally uniformly. Every contour deformation below changes only a bounded subarc about S_c and stays in this zero-free S-sector.

At the cusp the exact stationary equation and derivatives give

    Phi_S=Phi_SS=Phi_SSS=0,
    Phi_SSSS=-2/L^4,
    beta:=L Phi_Sw=-L/(d w_c)<0.

Use the linear field eta=beta*(w-w_c). It equals (w/w_c-1)/c_d with c_d=-d/L, exactly the normalized field of the retained bounded-cusp theorem. Write eta=rho^3 omega, rho=|eta|^(1/3)>0, |omega|=1. When Im(w)>0, Im(eta)<0. Put S=S_c+L rho u and

    f_(rho,omega)(u)=rho^(-4)*[Phi(S_c+L rho u,w)-Phi(S_c,w)],
    k=m rho^4.                                               (C1)

On every fixed complex u-disc the family extends analytically to rho0 and converges with all fixed derivatives to

    f_(0,omega)(u)=-u^4/12+omega*u.                            (C2)

This follows by a convergent Taylor series, using the quartic jets and the linear field normalization. Parameters omega may range in a small complex neighborhood of any compact portion of the unit lower semicircle. No assertion that rho=|eta|^(1/3) is holomorphic is needed: it is only a parameter for uniform estimates. The EXACT periods introduced below are holomorphic functions of w.

## 1. Complete real tails after the blow-up

At w_c the positive real cycle has a unique global quartic maximum S_c. There are r_0,c,C>0 such that, for real x=S-S_c with |x|<=r_0,

    Phi(S,w_c)-Phi(S_c,w_c)<=-c*x^4.

The variation after subtracting the central value is bounded by

    |[Phi(S,w)-Phi(S_c,w)]-[Phi(S,w_c)-Phi(S_c,w_c)]|
             <=C*|w-w_c|*|x|.                               (C3)

Indeed Phi_w=1/(w+S^d), and this function minus its value at S_c is uniformly Lipschitz on the fixed local real interval and a small complex w-disc. Therefore on the original real local interval the scaled exponent obeys

    Re f_(rho,omega)(u)<=-c_1*u^4+C_1*|u|.                    (C4)

For real u its imaginary part is bounded by C_1 |u|, since the base phase at w_c is real and the perturbation in(C3) is O(rho^4 |u|). Subtracting any selected critical value adds only a uniformly bounded number in scaled coordinates. For sufficiently large fixed R, (C4) consequently gives a strict quartic loss on |u|>=R, including a bound |Im(f-f_i)|<=C|u|+C. A small fixed sector for k, |arg k|<epsilon_m, preserves this loss. This treats the entire interval R L rho<=|S-S_c|<=r_0.

Outside |S-S_c|<r_0 the original positive phase has a fixed strict gap. Its logarithmic left endpoint and exponential right tail are integrable, also with fixed powers of the phase. The perturbation of Phi and of its central value is uniformly O(|w-w_c|) on the ENTIRE positive ray because w+S^d stays away from0 for Re(w)>0. For sufficiently small rho0, the original nonlocal gap remains fixed. For complex m with |arg m|<epsilon_m the imaginary perturbation is still bounded, so these whole tails are exponentially small relative to every selected packet. Their relative size after the Gaussian prefactor is at most C sqrt(|m|) rho exp(-c|m|), which is O(1/|k|) when |k|>=K and rho<=rho0. There is no hidden factor exp(C m |w-w_c|) left uncontrolled: the fixed base gap is preserved before multiplying by m.

## 2. Outer-angle localization by a single horizontal packet

Fix a small angular epsilon>0. Suppose first that Re(omega)>=sin(epsilon), Im(omega)<=0. The dominant saddle of(C2) is

    u_R=3^(1/3)*exp(i arg(omega)/3),

whose real part is positive and whose argument lies in[-pi/6,0]. Put y_R=Im(u_R). Along u=x+i y_R,

    Re f_0= -x^4/12+y_R^2*x^2/2+Re(omega)*x
                         -y_R*Im(omega)-y_R^4/12.           (C5)

The x-dependent function V satisfies V(x)-V(-x)=2 Re(omega)x>0 for x>0. Its derivative -x^3/3+y_R^2 x+Re(omega) has exactly one positive root: it starts positive, increases to x=|y_R| and then decreases strictly to minus infinity. That root is Re(u_R). It is the unique global maximum, since every negative x has smaller value than its positive reflection. Moreover

    -Re f_0''(u_R)=Re(u_R^2)>=3^(2/3)/2>0.

On the compact angle set in question this maximum has uniform quadratic curvature and a strict gap outside each small neighborhood. Reflection u->-conjugate(u), omega->-conjugate(omega) proves the corresponding statement for Re(omega)<=-sin(epsilon), with the left saddle.

For rho sufficiently small the analytic implicit-function theorem gives the corresponding actual simple saddle. On its horizontal line the single maximum, the quadratic bound and all compact strict gaps persist by C^k convergence in(C2). Choose R large enough that both vertical segments at Re(u)=+-R, from the real line to this horizontal line, have values below the saddle by a uniform gap. This is possible because the -R^4/12 term dominates the bounded-height corrections. The original real tails are already controlled by Section1. Cauchy deformation through this finite rectangle is legitimate on the narrow zero-free S-sector.

Expanding about the actual critical point on a symmetric real horizontal interval gives a nonzero complete Gaussian leading term. Its first cubic correction is odd and integrates to zero even for complex k and complex Gaussian curvature. The quartic and squared-cubic remainder has an integrable majorant C(|k| y^4+|k|^2 y^6) exp(-c |k| y^2), after shrinking the interval, so its integral is O(|k|^(-3/2)) against the |k|^(-1/2) leading term. The first relative error is O(1/k), uniformly in rho and in these outer angles. All other contour pieces have strict exponential loss. Thus the actual I_m(w) is NONZERO for all sufficiently large real k, uniformly in the outer-angle regions. This statement includes their possible Stokes directions; it imports no model zero-location estimate. Every actual zero with large k must lie in the smaller central cone

    |arg(eta)+pi/2|<epsilon.                                  (C6)

Use a wider central cone of width2epsilon in the next sections, leaving a margin for radial transport.

## 3. Exact analytic split and two uniform packets in the central cone

At omega=-i the two lower saddles of(C2) are

    u_R=r*(sqrt(3)/2-i/2), u_L=r*(-sqrt(3)/2-i/2), r=3^(1/3).

Their common horizontal line has two strict real-part maxima, and the midpoint -i r/2 lies below both by9r/64. Each halfline has exactly one of the maxima. In fact, writing c=r/2 and x_i=+-sqrt(3)c,

    Re[f_0(u_i)-f_0(x-i c)]=(x^2-3c^2)^2/12,
    Im[f_0(x-i c)-f_0(u_i)]=(c/3)*(x-x_i)^2*(x+2x_i).

On the corresponding halfline x*x_i>=0 the absolute imaginary difference is at most5 times the real loss. Consequently a sufficiently small FIXED sector for k preserves the local quadratic loss and every complete model halfline gap. The midpoint and the outer cuts have strict loss as well.

For rho small and omega in the wider central cone, the two actual simple saddles and their midpoint remain in fixed small neighborhoods of these model points in u-coordinates. Let S_L(w),S_R(w) be their physical values and let

    S_M(w)=(S_L(w)+S_R(w))/2.

They are holomorphic on the punctured central w-cone, by their Puiseux branches or by the retained distinguished inverse-branch theorem. Define EXACT endpoint periods P_L(m,w) from0 to S_M(w), and P_R(m,w) from S_M(w) to infinity, using the positive ray outside a fixed small neighborhood of S_c and finite local paths within the narrow zero-free sector. Equivalently split the positive integral at a fixed real point near S_c and add/subtract the same finite path to S_M. These definitions give

    I_m(w)=P_L(m,w)+P_R(m,w)                                 (C7)

for all m in the complex sector, and both functions are holomorphic in(m,w). The integral is path-independent in this chosen logarithm domain; no untwisted homology for general m is used.

For estimates deform the left and right paths to their respective actual saddle horizontals and connect them to S_M by short central paths. The perturbed midpoint gaps put each central connector below BOTH peaks. Outer vertical cuts and all real tails are controlled by Section1. All these estimates persist for a small fixed complex-m sector: on the compact local contour the ratio of imaginary phase change to the strict real loss is uniformly bounded (near the saddle both vanish quadratically), and(C4) controls the unbounded local real tails. The far original tails retain their fixed gap. This establishes uniform Gaussian expansions for the EXACT functions(C7):

    P_i(m,w)=exp(m Phi_i(w))*sqrt(2pi/m)*c_i(w)*(1+R_i(m,w)),
    c_i(w)=1/sqrt(-Phi_SS(S_i(w),w)),
    |R_i(m,w)|<=C/|k|, i=L,R.                               (C8)

Choose square roots with the orientations of the actual left and right contours. These are nonzero holomorphic factors. The expansions hold uniformly for |arg m|<epsilon_m, small rho, the wider central cone and |k|>=K. A direct Taylor proof uses a symmetric horizontal interval: the Gaussian term is explicit; the odd cubic term cancels, the remaining local error is O(k^(-3/2)) against the k^(-1/2) main term, and all remaining pieces have exponential loss. All constants are fixed before m or rho is varied.

For every real m and w in a slightly smaller cone, the SAME uniform bound holds on a genuine product of a complex-m disc of radius c m and a complex-w disc of radius c |eta|=c rho^3. Shrink c using the angular and sector margins. Cauchy's estimate applied to the ACTUAL holomorphic functions R_i defined by(C8), not to a formal O-symbol, yields

    |partial_m R_i|<=C/(m k),
    |partial_w R_i|<=C/(k rho^3).                            (C9)

The use of |eta| in these estimates is only a size description of the discs. It does not enter the holomorphic definitions of P_i,R_i.

## 4. The Puiseux amplitude ratio and exact zero equation

Choose the analytic cube root t=eta^(1/3) in the central cone. The implicit-function theorem applied after S=S_c+L t u gives

    S_i=S_c+t a_i(t),
    Phi_SS(S_i,w)=t^2 q_i(t),

with a_i,q_i holomorphic and q_i(0) nonzero. The two leading cubic roots differ by a fixed cube root of unity, so |q_L(0)|=|q_R(0)|. Hence c_R/c_L is a nonzero analytic unit in t, with a limit of absolute value1. Choose its holomorphic logarithm b_0. Then

    |Re b_0(w)|<=C rho, |partial_w b_0(w)|<=C/rho^2.           (C10)

The derivative bound follows from dt/dw=beta/(3t^2). No separate claim that c_i remains bounded as rho tends0 is made; each individual c_i is of order rho^(-1), and its common singular factor cancels only in the ratio.

Increase K so that |R_i|<1/4 and define the EXACT holomorphic functions

    b(m,w)=b_0(w)+log(1+R_R(m,w))-log(1+R_L(m,w)),
    D(w)=Phi_R(w)-Phi_L(w),
    E(m,w)=m D(w)+b(m,w).                                    (C11)

For real m, exp(E)=P_R/P_L. At every actual zero of I_m the two nonzero periods have ratio-1, so

    Re E(m,w)=0.                                            (C12)

No zero index or imaginary-logarithm branch choice is required. From(C9)--(C10),

    |b_m|<=C/(m k),
    |b_w|<=C[rho^(-2)+1/(k rho^3)],
    |Re b|<=C[rho+1/k].                                     (C13)

The same Puiseux expansions give constants c,C>0 such that throughout the wider cone

    Re(S_R-S_L)>=c rho, |S_R-S_L|<=C rho, |D|<=C rho^4.       (C14)

The first inequality follows from the positive real leading difference at eta=-i|eta| and then a fixed angular margin. Before applying Cauchy, enlarge the final comparison threshold K by a fixed factor relative to the packet-existence threshold, so that its product discs remain in the latter domain. At the end choose the claimed physical neighborhood strictly inside the radius on which all packet estimates were constructed, reserving a fixed collar for the short actual radial path. At a stationary point, Phi_S=0 implies S^d/(w+S^d)=(S-a)/d, hence w Phi_w=(a+d-S)/d. Subtraction and stationarity, which removes the root derivatives, give the exact envelope identity

    d*w*D_w=-(S_R-S_L).                                      (C15)

## 5. Compare the actual consecutive orders

Suppose x is a common root of the actual consecutive moment polynomials at orders N and n=N-1. Set w=x/N^d. The exact actual moment normalization gives I_N(w)=I_n((N/n)^d w)=0. It uses each order's own Gamma shape; no fixed-shape replacement occurs.

First suppose |k_N|=N |eta(w)|^(4/3)<=2K. The normalized field of the retained uniform bounded-cusp theorem has

    |h_N|=N^(3/4)|eta|=k_N^(3/4)<=(2K)^(3/4).

The fixed-cusp bounded-window theorem applies with the finite absolute window L_s=(2K)^(3/4)/|beta|; alternatively its stronger uniform-H version applies directly. It excludes the pair for all sufficiently large N. This is exactly a fixed window at this step. No growing-field statement is imported.

It remains that k_N>=2K. Section2 excludes the outer angles, so w lies in the smaller central cone. Define

    m=N-theta,
    w_theta=[N/(N-theta)]^d*w, 0<=theta<=1.                   (C16)

Since d is fixed and w stays near the fixed positive w_c, |eta(w_theta)-eta(w)|<=C/N. Its relative size is O(1/(N rho^3))=O(rho/k_N). Choose rho0 sufficiently small, K sufficiently large and N>=2. The full path stays in the wider central cone, rho remains comparable to its initial value, and k_theta>=K by the2K overlap. All estimates(C8)--(C15) apply on the entire real-m path.

At its upper endpoint(C12) and(C13) imply

    |Re D(w)|<=C[rho/N+1/(N k_N)].

By(C15), along(C16),

    dD(w_theta)/dtheta=-(S_R-S_L)/m,

so |Re D(w_theta)|<=C[rho/m+1/(m k_theta)] uniformly on the path. Differentiate the EXACT E, using(C15):

    dE(m,w_theta)/dtheta
      =-D-(S_R-S_L)-b_m+(d*w_theta/m)*b_w.                   (C17)

The main real part is at most-c rho. Its errors are bounded, using m=k/rho^4, by

    C[rho/m+1/(m k)+1/(m rho^2)+1/(m k rho^3)]
        <=C[rho^5/k+rho^4/k^2+rho^2/k+rho/k^2].              (C18)

Constants absorb the fixed d and w_c. First choose K large and then rho0 small so that(C18) is less than c rho/2 throughout the comparable path. Therefore

    d Re E/dtheta <=-c rho/2<0.

But both actual endpoint zeros require Re E=0 by(C12). This contradiction excludes every common root in the punctured small physical cusp neighborhood with k_N>=2K. Together with the retained bounded-field core theorem it proves:

**Theorem FixedCuspNeighborhoodConsecutiveExclusion.** For every fixed integer cusp there exist r_s>0 and an unevaluated N_s such that for all N>=N_s the actual consecutive moment polynomials have no common complex root x with |x/N^d-w_c|<r_s.

The lower half-plane follows by real coefficients and conjugation. On the positive real line the actual moment polynomial has positive coefficients and therefore no root. The exact center also has no root. The constants here need not be uniform in s. This theorem covers ALL scaled fields within the fixed physical neighborhood; no hidden bounded-index condition remains.

## Scope and next deduction

The proof establishes a local-in-the-physical-plane statement, much stronger than a fixed scaled window. Combining it with a genuinely proved compact complement theorem across Stokes walls and real boundaries, and the all-order annulus for every fixed support, would establish eventual Strong Factorial on each fixed integer cusp. That global conclusion is not asserted by this note until its separate complement hypotheses are checked.

All analytic estimates in the proof are written contour arguments. Exact finite algebraic controls can support their jets, scale exponents, and derivative identities, but do not replace the complete-cycle, holomorphy, Cauchy, or uniformity proofs.
