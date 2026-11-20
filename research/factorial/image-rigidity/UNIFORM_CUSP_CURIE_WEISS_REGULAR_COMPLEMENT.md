# Uniform Curie--Weiss regular complement for integer Gamma cusps

This is a scratch proof of the uniform regular-complement step. It does not assert an effective cutoff or an all-order Strong Factorial theorem. The Curie--Weiss identification itself is classical and already recorded in GAMMA_CURIE_WEISS_CRITICAL_BRIDGE.md; the issue here is uniform complete original contours when both the cusp and moment order vary. No approximation requires m/s to be small.

## 1. Actual normalization, compactification, and root annulus

For an integer s>=1 set

    a=s^2, d=2s+1, C=s(s+1), L=sqrt(C),
    kappa=(s+1)/s, w_c=C^d*kappa, q=w/w_c,
    t=(S-C)/L, p=1+t/L.

The original period at integer order m is, up to the positive nonzero factor L*[C^(a+d)*exp(-C)]^m, exactly

    J_(s,m)(q)=integral_(-L)^infinity Psi_s(t,q)^m dt,
    Psi_s(t,q)=p^a*exp(-L*t)*(kappa*q+p^d).             (G1)

For every integer s,m the finite-contour integrand is entire in t and polynomial in q. In particular, Cauchy deformations of a finite middle segment are legitimate even if that segment crosses a zero of Psi_s. No global logarithm or division by q+1 is used.

On each fixed complex t-disc, Psi_s and every fixed derivative converge, locally uniformly also in q, to

    Psi_0(t,q)=exp(-t-t^2/2)*(q+exp(2t))
             =exp(-t^2/2)*(q*exp(-t)+exp(t)).          (G2)

More precisely the family is analytic in sigma=1/s at sigma=0 on each such disc. Put kappa=1+sigma. The apparent singularity cancels:

    a*Log(p)-L*t
       =-t/sqrt(kappa)-t^2/(2*kappa)
        +sum_(j>=3) (-1)^(j+1)*sigma^(j-2)*t^j/
                                  [j*kappa^(j/2)],

and d*Log(p)=(2/sigma+1)*Log(1+sigma*t/sqrt(kappa)) is analytic at0 with limit2t. The bracket in (G1) need not be nonzero; (G2) is an identity and convergence of the entire functions. Local logarithms will be chosen only near nonzero critical factors and on contour pieces where they are needed.

The already proved sharp Gamma annulus gives, at every actual root of order m>=2,

    1/16 < |q| < 8.                                    (G3)

Indeed a^d/w_c=(s/(s+1))^(2s+2)>=1/16 and
(a+d)^d/w_c=((s+1)/s)^(2s)<exp(2)<8. Thus every pair argument lies in a single fixed compact annulus.

## 2. Uniform complete real tails

Let r range in any compact positive interval [r_-,r_+]. On the original real t domain p>0, Psi_s(t,r)>0 and

    H=p^d/(kappa*r+p^d), 0<H<1,
    partial_t log Psi_s(t,r)=(-s-L*t+d*H)/(L*p).         (G4)

Consequently

    t>=4:  partial_t log Psi_s(t,r)<=-t/(2p)<=-2/5,
    t<=-2: partial_t log Psi_s(t,r)>=(-t-1)/p>=1.        (G5)

For the first inequality use (s+1)/L=sqrt(kappa)<=sqrt(2) and p<=1+t. For the second use s/L<=1 and p<=1. These bounds apply throughout the original real tail, including its finite left endpoint t=-L; that endpoint causes no missing interval.

For complex q,

    |Psi_s(t,q)|<=Psi_s(t,|q|).                         (G6)

Given a fixed compact q-set away from0, choose a large fixed R>=4. For all sufficiently large s, L>R and (G5) implies

    integral_R^infinity |Psi_s(t,q)|^m dt
          <=(5/(2m))*Psi_s(R,|q|)^m,
    integral_(-L)^(-R) |Psi_s(t,q)|^m dt
          <=(1/m)*Psi_s(-R,|q|)^m.                     (G7)

On each fixed horizontal-height interval |Im t|<=Y,

    |Psi_0(+-R+i y,q)| <=
             exp[-R^2/2+Y^2/2]*(r_+*exp(R)+exp(R)).

This tends to0 as R tends to infinity. The same holds for the real endpoint values in (G7). First choose R, then s sufficiently large for convergence on the fixed rectangle. Thus both complete original tails and the two vertical connectors at +-R can be placed below ANY fixed positive critical-modulus threshold. Their relative contributions are O(exp(-c*m))/m, uniformly in all sufficiently large s. This is why there is no need to control the complex descent flow all the way to the receding endpoint -L.

Only real integer orders are needed in Sections3--7. In particular, these statements do not claim a global untwisted cycle identity for arbitrary complex m.

## 3. The two limiting physical saddles

At a nonzero critical factor of Psi_0, the stationary map is

    q=Q(t):=exp(2t)*(1-t)/(1+t),
    Q'(t)=-2*t^2*exp(2t)/(1+t)^2.                       (G8)

On the open physical strip

    -pi/2 < Im(t) < 0,

every q in the upper half-plane has exactly two simple inverse points, denoted t_L(q),t_R(q), with

    Re t_L<0<Re t_R.                                   (G9)

Here are the count and ordering details. The top boundary maps to the real line, apart from its boundary pole t=-1. On the bottom boundary exp(2t) is negative real and Im[(1-t)/(1+t)]>0, so Q maps strictly into the lower half-plane. The only finite zero of Q' is the boundary point t=0, with value1. For target q in a compact subset of the upper half-plane no preimage can escape horizontally: Q(t) tends to0 as Re t tends to minus infinity and has modulus tending to infinity as Re t tends to plus infinity. No preimage can approach the boundary pole at a bounded target value. Hence the inverse count is finite and constant in the upper half-plane.

At q=iM, M large, there is one inverse near the simple pole -1, namely -1+2*exp(-2)/q+O(q^-2), and one with large positive real part, t=(Log(-q))/2+O(1/Log M), taking Im Log(-q) in(-pi,0). There are no others: after removing a small pole disc, large targets force Re t large, where the logarithmic equation
t=Log(-q)/2+Log((1+t)/(t-1))/2 is a contraction on its forced O(1/Log M) neighborhood. Thus the count is2. Since the upper target half-plane is simply connected and has no branch value, both inverse branches are holomorphic.

For t=i y, 0<y<pi/2, the continuous argument of Q(i y) is

    2y-2*atan(y), strictly between0 andpi.

Therefore Q(-i y) lies in the lower half-plane, so neither upper-target inverse crosses Re t=0. The large-target ordering gives (G9) everywhere. On the negative target axis the two branches extend to the two distinct real roots t_L<-1<t_R, with t_R>1. On the positive axis other than1, the physically dominant real critical point is simple. All these assertions also follow directly from the strictly decreasing real map Q on(-1,1) and its two exterior negative-value intervals.

The critical factor at a stationary point is

    lambda_0(t)=2*exp(t-t^2/2)/(1+t),                     (G10)

so it is finite and nonzero on any compact target set under consideration. Analytic inverse continuation and (G8) keep the two saddles in a fixed t-rectangle on compact sets avoiding the cusp target q=1.

## 4. Exact horizontal maximum lemma

Write q=exp(2h+i theta), h real. On a horizontal line t=x+i y, put X=x-h and C_*=cos(2y-theta). Apart from a constant independent of X,

    log|Psi_0(x+i y,q)|
        =-X^2/2-h*X+(1/2)*log(cosh(2X)+C_*).            (G11)

At a zero of the bracket the modulus is0; all maximum statements below are made at nonzero factors.

For h>0, reflection favors X<0:

    V(-z)-V(z)=2h*z>0, z>0.

On that favored side, any stationary maximum has z=-X>h, and its equation is

    C_*=T_h(z):=sinh(2z)/(z-h)-cosh(2z).                 (G12)

The function T_h decreases strictly from plus infinity to minus infinity on z>h. To see this, put v=z-h>0. The numerator of -T_h' is

    (1+2v^2)*sinh(2z)-2v*cosh(2z)>0.

After division by cosh(2z), monotonicity of tanh reduces this to

    E(v)=(1+2v^2)*sinh(2v)-2v*cosh(2v)>0,
    E(0)=0, E'(v)=4v^2*cosh(2v)>0.                      (G13)

Thus each horizontal line has exactly one favored global maximum. It is simple: at a root of (G12),

    V''=(z-h)*T_h'(z)/(cosh(2z)+C_*)<0.                 (G14)

In particular the horizontal through t_L(q), whose real part is negative, has that critical point as its unique global maximum. For h<0 reflection gives the same statement with t_R.

For h=0, C_*<1, T_0 decreases strictly from1 to minus infinity on z>0, by (G13). Therefore the horizontal has exactly two symmetric simple global maxima at +-z; the middle is strictly lower (and can be a zero when C_*=-1).

If |q|=1, q is upper nonreal, the exact involution

    Q(-conjugate(t))=1/conjugate(Q(t)),
    Psi_0(-conjugate(t),q)=q*conjugate(Psi_0(t,q))       (G15)

interchanges t_L and t_R. They lie on the SAME horizontal line. By (G9) they are distinct and by the preceding calculation they are the two simple global maxima of the original whole limiting horizontal. The negative endpoint q=-1 has the same statement, with both critical points real and C_*=-1. On every closed unit-circle arc avoiding q=1, the saddle separation, the negative real Gaussian curvature, the middle loss and every outer loss have uniform positive margins.

This is a direct complete-contour fact about the limiting phase, not a transfer of a finite-degree circle law to the Gamma family.


## 5. Complete actual one-packet localization off the unit circle

Fix any nonzero q_* with |q_*|!=1. If q_* is upper nonreal, use the unique dominant limiting horizontal from Section4. For a positive real q_* use its real simple maximum; for a negative real q_* use the favored maximum of the real horizontal with C_*=-1. The lower half-plane follows by conjugation.

Choose R large enough for Section2. On the finite limiting horizontal [-R,R]+iy, its dominant maximum has nonzero factor, strictly negative real quadratic curvature, and a strict gap on the rest of the contour outside a small saddle disc. The two vertical pieces have strict loss, as do the whole original real tails by (G7). All inequalities persist in a sufficiently small complex q-neighborhood and for all sufficiently large s: use convergence only on this FIXED finite rectangle, and keep the ACTUAL saddle, factor and curvature in the leading term.

A symmetric horizontal interval through the actual saddle gives

    J_(s,m)(q)=lambda_s(q)^m*m^(-1/2)*
                       c_s(q)*(1+O(1/m)),                  (G16)

uniformly in that neighborhood, all sufficiently large s and all integer m. The analytic nonzero c_s is sqrt(2pi)/sqrt(-partial_t^2 Log Psi_s) with its contour orientation. The first cubic correction is odd on the symmetric interval; the fourth and squared-third corrections have Gaussian-integrable majorants, giving the stated relative O(1/m). All remaining pieces have strict exponential loss. This derivation uses actual Psi_s in lambda_s and c_s. It makes no estimate of exp(m*Phi_s) by exp(m*Phi_0), so m/s is unrestricted.

It follows by a finite cover of any compact set separated from the unit circle that the actual J_(s,m) is zero-free there for all sufficiently large s,m. Equivalently:

**Uniform asymptotic circle capture.** If s_j tends to infinity, m_j tends to infinity, and q_j are actual normalized roots, then every subsequential limit of q_j has modulus1.

The annulus (G3) supplies compactness. This is not an exact circle law for any finite Gamma polynomial; those already fail a common-modulus test in degree3.

## 6. Two complete actual packets near the limiting unit circle

Fix q_* on the unit circle with q_*!=1. Treat upper and lower conjugate arcs; the q_*=-1 neighborhood is allowed. Section4 gives two strict simple maxima on a common limiting horizontal, with a strict middle gap. Choose R and a small q-disc U about q_* so that both critical factors are nonzero and both saddles extend analytically on U.

For all sufficiently large s the two actual saddles t_(s,L),t_(s,R) are jointly analytic perturbations, with

    Re(t_(s,R)-t_(s,L))>=c>0, |t_(s,i)|<=C.               (G17)

The normalized phases, factors, Gaussian leading constants and each fixed q-derivative are uniformly bounded on a slightly larger disc; critical factors and Gaussian curvatures have uniform positive lower bounds.

For each integer m split the EXACT original period at the analytic midpoint M_s(q)=(t_(s,L)+t_(s,R))/2. Outside [-R,R] keep the original real tails. Inside use finite paths from -R to M_s and M_s to R. Since the integer integrand is entire, this gives an exact sum J=P_L+P_R, independent of finite path choices. For estimates deform each half to its own saddle horizontal with short connectors to M_s. At the limiting parameter every middle connector lies below BOTH peaks; the same remains true after shrinking U and taking s large. Vertical outer connectors and the whole real tails have the uniform loss established earlier.

Consequently

    P_(i;s,m)(q)=lambda_(s,i)(q)^m*m^(-1/2)*
                           c_(s,i)(q)*(1+R_(i;s,m)(q)),
    |R_(i;s,m)|<=C/m, i=L,R.                              (G18)

All constants are uniform in large s, integer m and q in U. At q_*=-1 the midpoint can coincide with a zero of Psi; that is harmless for the integer entire integrand and lies in the strict-gap region. Logarithms are only chosen near the two nonzero critical factors. No complex-order identity for a zero-crossing middle path is required.

Choose local analytic logs and write

    D_s=Log lambda_(s,R)-Log lambda_(s,L),
    b_(s,0)=Log c_(s,R)-Log c_(s,L).

At an actual zero, (G18) gives the value estimate

    Re[m D_s(q)+b_(s,0)(q)]=O(1/m).                       (G19)

The derivatives of the ACTUAL leading amplitude satisfy |partial_q b_(s,0)|<=C. This follows from local analytic curvature and nonvanishing, uniformly in sigma=1/s. No derivative of R_i or of an unspecified O-symbol is used.

The actual original stationary identity w*Phi_w=(a+d-S)/d yields the exact normalized envelope equation

    q*D_s'(q)=-(L/d)*(t_(s,R)-t_(s,L)).                   (G20)

The scalar used to normalize Psi_s is independent of q and cancels from the action difference.

## 7. Uniform actual consecutive exclusion on the regular complement

For every r>0, there is an unevaluated N_r such that for EVERY integer cusp s and every N>=N_r, the actual consecutive pair has no common root with

    |q_N-1|>=r, q_N=x/(N^d*w_c).                         (G21)

Here no individual upper or lower moment is asserted zero-free near the unit circle. The conclusion concerns their actual consecutive common roots.

Suppose otherwise. There is a sequence N_j tending to infinity with common roots violating (G21). If s_j has a bounded subsequence, the already proved fixed-support eventual theorem excludes it. Thus after passage to a subsequence s_j tends to infinity. Let n_j=N_j-1 and

    q_(n_j)=(N_j/n_j)^d*q_(N_j)=lambda_j*q_(N_j),
    lambda_j>1.

Both endpoint roots lie in (G3), and Section5 gives |q_N|,|q_n| tending to1. Therefore lambda_j tends to1, both arguments converge to the SAME unit-circle point q_* after passage to a subsequence, and q_*!=1 by (G21). Also d/N tends to0 because log(lambda_j)=d*log(N/n).

Choose the uniform two-packet disc U from Section6. The ENTIRE actual adjacent path

    m=N-theta, q_theta=(N/m)^d*q_N, 0<=theta<=1,

lies in U for all sufficiently large j. Let Delta_s=t_(s,R)-t_(s,L) and

    E_0(theta)=m*D_s(q_theta)+b_(s,0)(q_theta).

Use (G20), not a frozen limiting action, to differentiate:

    E_0'=-D_s-L*Delta_s+(d*q_theta/m)*b_(s,0)'(q_theta).  (G22)

At the upper root, (G19) implies |Re D_s(q_N)|<=C/N. Along the path, dD_s/dtheta=-L*Delta_s/m, so

    |Re D_s(q_theta)|<=C*(1+L)/N.

Since d/L=(2s+1)/sqrt(s(s+1)) is uniformly bounded, (G17), (G22) and the uniform derivative bound on b_(s,0) give, for all large N,

    d Re E_0/dtheta <=-c*L+C*(1+L+d)/N <=-c*L/2.        (G23)

Both endpoint roots, however, require Re E_0=O(1/N) by (G19). Integrating (G23) contradicts that, since L>=sqrt(2). This proves (G21).

This argument handles regular Stokes directions without a new Picard--Lefschetz import: the explicit limiting horizontal is either a strict one-peak contour or a strict two-peak contour, and all actual finite central contours are its uniform perturbations.

## 8. Relation to the remaining local problem

Section7 is the GLOBAL regular-complement input. A separate uniform physical cusp-neighborhood theorem would assert that there exist r_*>0,N_* independent of s with no actual common root in |q_N-1|<r_* for N>=N_*. Combining that theorem with (G21) gives ONE unevaluated N_0 valid for ALL integer cusps and all coefficients/characteristic-zero fields, by the usual rational Bezout equivalence.

For that local proof the central field is exactly

    eta=-(L/d)*(q-1), rho=|eta|^(1/3), k=m*rho^4.

Near q=1 the compact normalized phase is analytic through sigma=0. Its fourth t-derivative is -2, and its raw mixed second-t/first-q derivative vanishes at the cusp. Thus the uniform quartic packet argument uses a common local analytic domain and the real tails from Section2 or the canonical402 concentration proof.

A useful exact path observation is as follows. If BOTH actual endpoint zeros are in the central lower eta-cone, then

    eta_theta=lambda_theta*eta_N-(L/d)*(lambda_theta-1),
    Im eta_theta=lambda_theta*Im eta_N,

and Re eta_theta/Im eta_theta is monotone between its endpoint values. Hence the whole radial path remains in that cone. An upper central large k already forces lower k comparable from its imaginary part before applying the lower cone localization. One must first ensure that the lower endpoint is inside the local physical disc; this follows either from a uniform positive-real one-packet catch or, in a contradiction sequence, from Section5, which forces q_n to the same limiting unit-circle point1. The observation alone is not a lower-disc capture theorem.

If this local input is supplied, the set of bad pairs (s,N) is FINITE, not merely finite at every fixed s. Indeed N<N_0 at a bad pair; the proved sharp soft criterion excludes every pair with s^2>=(2s+1)(N-2), and its reverse inequality bounds s for each of the finitely many remaining orders. No assertion that this finite set is empty, or an effective bound on it, follows.

## 9. Antecedents and scope

The exact limiting polynomial is, up to a nonzero Gaussian scalar,

    sum_(k=0)^m binom(m,k)*q^(m-k)*
                              exp[(2k-m)^2/(2m)],

since direct Gaussian integration of (G2)^m gives that expression. This is the finite critical Curie--Weiss model at beta=1, hence has the classical finite Lee--Yang circle law. The same identification with a unitary Hermite polynomial, and the genuine O(m/d) obstruction to coefficient-relative transfer, were already recorded in GAMMA_CURIE_WEISS_CRITICAL_BRIDGE.md; see Kabluchko arXiv2203.05533v3, equations(1.5),(1.7)--(1.9), Lemma2.1 and Theorem2.18. Canonical QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md Section3 and UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md already use the critical CW/quartic limit.

None of these primary circle theorems is a premise of Sections4--7. The direct horizontal maximum calculation supplies the required limiting contour geometry. The actual leading exponent is always Phi_s at its ACTUAL stationary point. The uniform geometric convergence is used only for strict gap, curvature and finite-contour stability; no assumption m/s=O(1) or m/s tending to0 is made.

This note is a research proof, not an assertion of independent replay, external review or formal verification. Finite algebra checks support identities but cannot replace the whole-contour or uniformity arguments.

