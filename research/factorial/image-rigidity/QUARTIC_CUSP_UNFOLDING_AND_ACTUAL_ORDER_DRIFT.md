# Quartic cusps of the actual ordinary-binomial phase

This research note proves the local analytic and algebraic statements below. It does not prove consecutive-moment coprimality near a cusp. The gap d is an integer at least two. Local parameters may be complex; statements about the original positive integral specify real parameters.

## 1. The actual phase and its cusp

Put delta=d/sqrt(a), p=1+delta*v/d, and

    h(v)=(p^d-1)/delta,
    U(v)=d^2/delta^2*(log(p)-delta*v/d),
    F(v,z,delta)=U(v)+log(z+h(v)).

Branches are positive real at the cusp and locally holomorphically continued. The substitution S=a+sqrt(a)*v, w=a^d*(-1+delta*z) takes the original m-th moment integral to integral exp(m*F)dv times the nonzero scalar sqrt(a)*[a^(a+d)*delta*exp(-a)]^m.

Write

    R=(d+1)/(d-1), delta_c=R+1=2d/(d-1),
    a_c=(d-1)^2/4, v_c=1,
    z_c=(R^(d+1)+1)/delta_c.

At the cusp p=R and z_c+h(1)=R^d. The physical stationary point is S_c=(d^2-1)/4, and w_c/a_c^d=R^(d+1)>0. Integer cusps are d=2s+1,a=s^2,S_c=s(s+1), s>=1.

**Lemma QuarticCuspJets.** At this point,

    F_v=F_vv=F_vvv=0,
    F_vvvv=-2/R^2,
    F_vvvvv=24/[(d-1)*R^3].                     (Q1)

The parameter-jet matrix for (z,delta) is

    F_vz=-R^(-d-1),             F_vvz=0,
    F_vdelta=(R+2-R^(-d))/[R*(R+1)^2],
    F_vvdelta=(d-1)^3/[d*(d+1)^2].              (Q2)

Its determinant is

    -R^(-d-1)*(d-1)^3/[d*(d+1)^2].             (Q3)

Uniformly for d>=2,

    -F_vvvv>=2/9,       |det(Q2)|>=1/486.       (Q4)

Proof. The exact stationary numerator is

    P=z*v+(v-delta)*h-1,
    F_v=-P/[p*(z+h)].                          (Q5)

At the cusp h'=R^(d-1),h''=2R^(d-2),

    h'''=4(d-2)/(d-1)*R^(d-3),
    h''''=8(d-2)(d-3)/(d-1)^2*R^(d-4).

Thus P=P'=P''=0, P'''=2R^(d-1), P''''=8(d-2)/(d-1)*R^(d-2). The logarithmic derivative of p*(z+h) is exactly one there. Differentiating (Q5) gives (Q1). The z derivatives follow from F_vz=-h'/(z+h)^2. At fixed v,z,

    h_delta=[delta*v*p^(d-1)-p^d+1]/delta^2,
    (h')_delta=(d-1)*v*p^(d-2)/d.

In F_vvdelta the h_delta terms cancel because h''/(z+h)=2/R^2 and (h')^2/(z+h)^2=1/R^2. Simplification gives (Q2)-(Q3). Finally (d-1)^3/[d*(d+1)^2]>=1/18, while R^(d+1)<=27: differentiating (x+2)log(1+2/x), x=d-1>=1, gives log(1+2/x)-2/x<0. This proves (Q4).

## 2. A self-contained uniform local normal form

**Lemma UniformLocalQuarticNormalForm.** There are r,eta>0 independent of d>=2 such that for |z-z_c|+|delta-delta_c|<eta there are holomorphic C,H,T and an invertible local coordinate v=V(t,z,delta;d), |t|<r, satisfying

    F(V(t),z,delta)=C+H*t+T*t^2/2-t^4/12.      (Q6)

The functions are real on real parameters, H=T=0 at the cusp, and V_t(0)>0 there. These are existence radii, not numerical remainder bounds.

Proof. Compactify d by nu=1/d in [0,1/2]. The expressions

    p=1+nu*delta*v,
    h=[exp(log(p)/nu)-1]/delta,
    U=[log(p)-nu*delta*v]/(nu^2*delta^2)

extend holomorphically through nu=0, where h=(exp(delta*v)-1)/delta and U=-v^2/2. The cusp centres delta_c=2/(1-nu),R=(1+nu)/(1-nu),z_c are jointly analytic near this compact segment, using the removable expression for R^(1/nu+1). On it p=R>0 and z_c+h=R^(1/nu)>0; hence one uniform neighbourhood supports the logarithms.

For x=v-1 the analytic quotient

    J_nu(x)=-12*(F(1+x,z_c,delta_c)-F_c)/x^4

has J_nu(0)=R^(-2)>=1/9. Compactness gives a common x-disc where it is nonzero. The positive fourth root defines t=x*J_nu(x)^(1/4), with a joint inverse V_0 on a common smaller disc. The base phase becomes exactly F_c-t^4/12.

Here is a parameter proof rather than an asymptotic-theorem import. On the disk algebra A_r seek V=V_0(t+u(t),nu) and C,H,T making (Q6) exact, with r smaller than the common inverse radius. Composition is analytic for small ||u||. The linearized equation at the base is the bounded isomorphism

    (u,C,H,T) |-> -t^3*u/3-C-H*t-T*t^2/2

from A_r x C^3 to A_r. Its inverse removes the first three Taylor coefficients and divides the remainder by t^3, bounded by Cauchy and the maximum principle. The analytic Banach implicit-function theorem gives the local solution. This derivative is independent of nu, and all higher derivatives are bounded on a compact family of slightly larger discs, so a common eta and common solution ball work. Uniqueness patches the solutions. Conjugation gives reality; shrinking the domains keeps V_t nonzero. This proves the lemma and supplies no actual-period remainder or global contour claim.

At nu=0, with x=v-1,

    F(1+x,z_c,2)-F_c=log(cosh(x))-x^2/2,
    F_c=3/2.                                 (Q7)


## 3. The canonical parameter dictionary

At the cusp the normalized base inverse has

    V(t)=1+sqrt(R)*t+3*t^2/[5(d-1)]+O(t^3),
    A=V_t(0)=sqrt(R), B=V_tt(0)=6/[5(d-1)].    (Q8)

Compare the fifth coefficient in (Q1) with F-F_c=-t^4/12 and invert to obtain (Q8). Parameter-dependent coordinate changes contribute first at degree three, since F_v vanishes to degree three. Thus the first parameter jets of (Q6) are

    H_z=-R^(-d-1/2),
    T_z=-6*R^(-d-1)/[5(d-1)],                 (Q9)

and H_delta=sqrt(R)*F_vdelta, T_delta=R*F_vvdelta+B*F_vdelta. In particular

    T_z/H_z = B/A = 6/[5*sqrt(d^2-1)].        (Q10)

The canonical parameter Jacobian is R^(3/2) times (Q3), so it too is uniformly nondegenerate. It is incorrect to infer T_z=0 from raw F_vvz=0: the nonzero V_tt term contributes. This matters for Lee-Yang statements requiring real T.

## 4. Actual adjacent-order transport

Let N>=2 be an integer, n=N-1, and

    I_m(w)=integral_0^infinity [S^a(w+S^d)]^m exp(-m*S)dS.

Substituting S=(n/N)*r gives exactly

    integral [S^a(w+S^d)]^n exp(-N*S)dS
      =(n/N)^((a+d)*n+1)*I_n((N/n)^d*w).       (Q11)

Hence the actual consecutive moments in common N scale compare at

    z_n=z+alpha*(z-1/delta),
    alpha=(N/n)^d-1>0.                        (Q12)

At the cusp the radial first jets are

    dH/dalpha=-sqrt(R)/delta_c,
    dT/dalpha=-3/(5*d).                       (Q13)

Use z_c-1/delta_c=R^(d+1)/delta_c and (Q9). For fixed d as N tends to infinity, at upper parameter exactly the cusp, the scaled variables h_m=m^(3/4)*H, tau_m=sqrt(m)*T obey

    h_n-h_N=-sqrt(d^2-1)/2*N^(-1/4)+O_d(N^(-5/4)),
    tau_n-tau_N=-3/5*N^(-1/2)+O_d(N^(-3/2)).   (Q14)

Only these asymptotics are fixed-d; (Q11)-(Q13) are exact and require no d/N smallness.

For bounded scaled field and fixed delta=delta_c, (Q10) gives

    tau_m=b*m^(-1/4)*h_m+O_d(m^(-1)),
    b=6/[5*sqrt(d^2-1)].                      (Q15)

Thus tau generally becomes complex on the actual complex-z curve. The real field drift is transverse to imaginary h-axis model zeros, but this does not yet provide an actual-period zero band.

## 5. The original real cycle at the exact cusp

For real S>0 and w=w_c>0, phi(S)=a_c*log(S)+log(w_c+S^d)-S has exactly one stationary point S_c, its global maximum. For a_c<S<D=a_c+d put

    W(S)=S^d*(D-S)/(S-a_c).

Exact differentiation yields

    W'(S)=-d*S^(d-1)*(S-S_c)^2/(S-a_c)^2.      (Q16)

Thus W decreases strictly from infinity to zero, except for its horizontal point at S_c. The sign of phi' is positive before S_c and negative after it, including S<=a_c and S>=D. The positive real cycle therefore passes through the quartic maximum. The exact positive w_c is not a moment zero. A complex deformation with a uniform relative quartic-vector remainder is a separate obligation.

## 6. Classical source scope and the unresolved bridge

The primary NIST DLMF Section36.12(i), equations36.12.1-36.12.10, was read. It treats oscillatory integrals exp(i*k*f), real-analytic phase and smooth amplitude. K=2 gives three coalescing critical points and a quartic normal form, matched by critical values; its approximation uses the canonical integral and both first two parameter derivatives. The displayed approximation is not imported here as a relative-error theorem for the actual Laplace cycle. Source: https://dlmf.nist.gov/36.12 .

The primary publisher abstract for F.Ursell, Integrals with a large parameter. Several nearly coincident saddle points, Proc.Camb.Phil.Soc.72(1972),49-65, explicitly starts after making the phase polynomial and addresses uniform remainder inequalities. The full theorem was not obtained, so no precise Ursell remainder estimate is invoked. Source: https://doi.org/10.1017/S0305004100050945 .

Levinson's1960 canonical-form announcement was located in the AMS primary catalogue, but its PDF returned403. No uninspected theorem is used; Section2 proves the needed local assertion. Source: https://doi.org/10.1090/S0002-9904-1960-10397-7 .

For the canonical real quartic

    Z4(h,tau)=integral_R exp(-u^4/12+tau*u^2/2+h*u)du,

exact identities are Z_hhh=3*tau*Z_h+3*h*Z and Z_tau=Z_hh/2. Arithmetic's independent Newman-primary audit gives imaginary simple h zeros for real tau. That is model information; (Q15) violates the real-temperature premise at nonreal z.

A full cusp SFC proof still requires:

1. The actual positive-cycle deformation, with nonlocal bounds uniform in d and the complex parameter window.
2. A relative remainder for the three-component basis Z4,Z4_h,Z4_hh. Analytic division by -t^3/3+T*t+H does not alone give a global bound.
3. A zero band for the complex temperature curve (Q15), including the analytic amplitude. Simple model zeros give local estimates per fixed zero; an all-index band needs quantitative derivative ratios and flow.
4. Comparison with exact (Q12)-(Q13), plus localization of every possible actual common zero. Local normal form alone does not give such localization.

No all-order ordinary-binomial assertion near a_c is claimed.

