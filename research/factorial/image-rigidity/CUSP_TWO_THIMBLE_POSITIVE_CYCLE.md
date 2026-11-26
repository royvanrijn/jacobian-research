# Two distinguished thimbles for the actual positive-cycle cusp integral

Fix integers s>=1,a=s^2,d=2s+1. Throughout Im(w)>0. Let

    g_w(S)=S^a*(w+S^d)*exp(-S),
    Phi=a*log(S)+log(w+S^d)-S, f=Re(Phi)=log|g_w|,
    C={S: -pi/d<arg(S)<0}, alpha=pi/d.

Define the actual moment polynomial A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k). Its own-order integral normalization is I_m(w)=(a*m)!*A_m(m^d*w)/m^((a+d)*m+1).

The logarithms need only be local. Their derivative, f, and g_w are single-valued. The actual ordinary integrand g_w(S)^m dS is entire when m is a positive integer. No untwisted relative-homology statement for arbitrary noninteger m is asserted.

The retained theorem DistinguishedCuspSaddleSectorBranches in CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md proves that C contains exactly two stationary points S_L,S_R, both simple, with Re(S_L)<s(s+1)<Re(S_R). Its hypotheses are exactly this integer cusp and Im(w)>0. There is exactly one zero of w+S^d in C,

    zeta=(-w)^(1/d), arg(zeta)=(arg(w)-pi)/d.

Write lambda_i=g_w(S_i). They are nonzero. A non-Stokes parameter in this note means lambda_R/lambda_L is not positive real. This sufficient condition excludes every trajectory between the two critical points; it is not asserted to classify which potential Stokes parameters actually have such a trajectory.

## 1. The actual inward descent and all three ends

**Lemma CuspWedgeInwardDescent.** On the closed wedge, away from the critical points and the ends0,zeta,infinity, normalized downward flow is the rational vector field

    dS/dt=-1/Phi'(S).                         (T1)

It satisfies dPhi/dt=-1 locally, hence

    f(S(t))=f(S(0))-t,
    g_w(S(t))=g_w(S(0))*exp(-t).              (T2)

It is strictly inward on both angular boundary rays.

Proof. In the Euclidean metric the complex form of grad(f) is conjugate(Phi'), so -grad(f)/|grad(f)|^2=-1/Phi'. To check the boundaries put

    B(S)=S*Phi'(S)=a-S+d*S^d/(w+S^d).

On the upper ray S=r>0,

    Im(B)=-d*r^d*Im(w)/|w+r^d|^2<0.

On the lower ray S=r*exp(-i*alpha),

    Im(B)=r*sin(alpha)+d*r^d*Im(w)/|w-r^d|^2>0.

Since partial_theta f=-Im(B), negative gradient points to smaller angles on the upper ray and larger angles on the lower ray. Multiplying it by a positive factor preserves strict inwardness. These are the actual finite-w formulas, not a limiting-saddle argument.

At the ends one has, with positive logarithmic multiplicities,

    f=a*log|S|+O(1) near0,
    f=log|S-zeta|+O(1) nearzeta,
    f=(a+d)*log|S|-Re(S)+O(|S|^(-d)) atinfinity. (T3)

Because cos(arg(S))>=cos(pi/d)>=1/2, the last expression tends uniformly to minus infinity in the closed wedge. Radial derivatives are positive near each finite zero and negative at large radius. Thus for a sufficiently negative number T, the sublevel {f<=T}, including the ends, has exactly three contractible components: a sector cap at0, an interior disk atzeta, and an outer sector cap atinfinity. To see there are no other components, first fix disjoint end neighbourhoods with these radial properties, then choose T below the minimum of f on the remaining compact set.

The closed wedge compactified at0 andinfinity is a closed disk X. Explicitly, if t=r/(1+r), the map

    (r,theta) |-> (cos(pi*t),sin(pi*t)*(1+2theta/alpha))

identifies it with a closed disk, collapsing r0 and rinfinity to its two boundary endpoints. The point zeta is an interior point. Set E={0,zeta,infinity}; f=-infinity there in the extended sense. This compactification and(T3) make every finite closed f-band away from critical points compact and disjoint from E.

## 2. No boundary topology change; precisely two one-handles

**Lemma CuspBandAndHandleDeformation.** If a finite interval [u,v] contains no critical value of f, then X_v={f<=v} strongly deformation retracts to X_u. At either simple interior stationary point, passage through its critical value attaches one one-handle, with local core the descending thimble.

Proof. On the compact band u<=f<=v the vector field(T1) is smooth and bounded and has no pole. Flow until time f(S)-u for a point with f(S)>u, and fix X_u. Equation(T2) gives the stopping time exactly. The trajectory cannot exit through either angular boundary because of strict inwardness; it cannot reach an end in a finite band by(T3). ODE continuation on the compact band supplies the whole stopping time. This gives a continuous strong deformation retraction, including points of the angular boundary. A tangential critical point of f restricted to a boundary ray therefore causes no topology change: the full gradient is nonzero there and the same inward deformation applies.

At a simple critical point p, Phi-Phi(p) has a double zero and a nonzero quadratic coefficient. Taking an analytic square root of its nonzero unit gives a local coordinate u=x+iy with

    Phi=Phi(p)-u^2/2,
    f=f(p)-(x^2-y^2)/2.                      (T4)

Before the level f(p), the local sublevel has two x-lobes; after it, the lobes are joined by a band with horizontal core. The explicit hyperbola inequalities x^2-y^2>=+-2epsilon identify the added piece with a one-handle. Outside a smaller critical neighbourhood the preceding band flow identifies the sublevels, so no additional piece is added globally. The analytic coordinate pulls the Euclidean metric back to a positive conformal multiple; descending trajectories through p are exactly the real u directions. If the two critical values of f coincide, choose disjoint neighbourhoods and make the two local attachments simultaneously; the argument is unchanged.

Start at the three contractible end caps and increase the level above the maximum of f on X. There are exactly two interior critical points, hence exactly two one-handles, and no boundary changes. The final X is connected and contractible. Each one-handle can reduce the number of connected components by at most one; therefore both must connect distinct current components. Neither can create a loop. The two descending handle cores consequently join the three end components as a tree.

## 3. Complete gradient thimbles and a relative basis

**Theorem CuspTwoThimbleTreeBasis.** For non-Stokes w, each critical point has a complete downward thimble Gamma_i consisting of its two descending arms. Every arm ends at one of0,zeta,infinity. The two thimbles join distinct pairs of ends and form a tree on E. Their oriented classes are a Z-basis of H_1(X,E;Z), which has rank2.

Proof. Begin each arm at a small nonzero real u in(T4), then use(T1). It cannot reach a side boundary in forward time, since the vector field is strictly inward. If continuation failed in a finite f-band, compactness would force a limit at a critical point. It cannot return to its initial critical point because f has strictly decreased. It cannot reach the other critical point because(T2) preserves arg(g_w), while the two critical factors have different arguments. All other limits in a finite band are regular and permit ODE extension. Thus every arm continues to f tending to minus infinity. Eventually it lies in one of the three end caps; nested deeper caps force convergence to that end. The handle cores of Section2 extended through these caps are the actual arms, so their tree incidence is as stated.

For completeness, the exact sequence of the pair gives

    H_1(X,E;Z) isomorphic to {(n0,nzeta,ninfinity) in Z^3: sum n=0},

because X is contractible. The boundary vectors of the two oriented edges of a three-vertex tree form an integral basis of this lattice. This proves the relative-basis assertion, not merely its rank.

A non-Stokes set is open and dense. The stationary envelope identity gives, on a local logarithm D=log(lambda_R/lambda_L),

    w*D_w=-(S_R-S_L)/d !=0.                  (T5)

The two factors have global holomorphic nonzero branches in the upper half-plane, so D may also be lifted there. The possible Stokes sets Im(D)=2pi*k have no open subset. Locally they are smooth real curves because D_w!=0. The tree and its oriented end labels are locally constant on non-Stokes flow chambers.


## 4. Absolute thimble integrability and endpoint caps

Near0, Phi'=a/S+O(1), so(T1) gives S_dot=-S/a+O(S^2). Nearzeta, Phi'=1/(S-zeta)+O(1), so (S-zeta)_dot=-(S-zeta)+O((S-zeta)^2). At infinity,

    Phi'=-1+(a+d)/S-d*w/[S*(w+S^d)] tends to -1. (T6)

For a sufficiently large radius, |Phi'+1|<=1/4, hence |S_dot|<=2 and |S_dot-1|<=1/3. Consequently radius is strictly increasing on an infinite arm eventually, because Re(S_dot*conjugate(S)/|S|)>=cos(alpha)-1/3>=1/6. Finite-end radii are strictly decreasing eventually by their displayed local equations. Thus each arm can be truncated by one small or large circular cut near its terminal end.

For a thimble arm, t=0 denotes the critical limiting endpoint; the ODE is used for t>0. Equation(T2) gives |g_w|=|lambda_i|*exp(-t). Its weighted length is integrable for every real m>0: near0, |S_dot|=O(exp(-t/a)); nearzeta, |S_dot|=O(exp(-t)); at infinity |S_dot|<=2. In particular every positive-integer thimble period

    P_i(m,w)=integral_(Gamma_i) g_w(S)^m dS    (T7)

is absolutely convergent.

The passage from relative topology to exact integrals needs an end argument. An integral over a finite sublevel is not zero merely because that sublevel is used in relative homology. Instead work with the actual ends E. On a radius-r cap at0 the integrand is O(r^(a*m)), so cap integrals are O(r^(a*m+1)); atzeta they are O(r^(m+1)). At infinity, for all sufficiently large R,

    |g_w(S)|<=2R^(a+d)*exp(-R/2)

uniformly on the closed sector arc, whose length is at most alpha*R. Its integral therefore tends to zero as R tends to infinity. Truncate any of our thimbles and endpoint paths, join their ends by these small circles or the large sector arc, and apply Cauchy's theorem to the finite closed curves. The integrand is entire, including at0 andzeta. Letting the caps shrink or expand proves that exact periods depend on their end-relative classes for these admissible curves. This justifies the use of H_1(X,E) without assuming a differential vanishes on a finite cap.

## 5. Exact original-cycle decomposition and activation scope

**Theorem ActualPositiveCycleTwoThimbleDecomposition.** For non-Stokes w and every positive integer m,

    I_m(w)=integral_0^infinity g_w(S)^m dS
           =n_L*P_L(m,w)+n_R*P_R(m,w),        (T8)

where n_L,n_R belong to {0,1,-1}, are independent of m and locally constant on each non-Stokes flow chamber. They are not both zero. The nonzero coefficients are exactly the oriented edges on the unique tree path from0 toinfinity. That path can have one edge or two.

Proof. The original positive ray is an admissible path on the upper boundary of X, with end boundary infinity-minus0. In a tree there is a unique integral edge chain with this boundary. It uses the path from0 toinfinity, each edge once with a sign, and all other edges with coefficient zero. The relative-basis theorem identifies this chain with the original relative class, and Section4 gives the exact integral identity. The tree and end labels are geometric, independent of m. Local constancy follows from continuous dependence of the finite flow segments on parameters, until they enter their stable end caps; non-Stokes separation excludes a critical-to-critical obstruction.

This theorem proves that only the two distinguished wedge saddles are needed. It does not prove both are active at every parameter. A rank-two homology group alone would not establish the coefficient statement, which here uses the actual two-handle tree and the original endpoint boundary.

The three fixed endpoint periods are useful concrete representatives:

    E_(0,zeta;m)(w)=integral_0^zeta g_w(S)^m dS,
    E_(zeta,infinity;m)(w)=integral_zeta^infinity g_w(S)^m dS,
    I_m=E_(0,zeta;m)+E_(zeta,infinity;m).

Each thimble joins two distinct ends, so its exact period is, up to orientation, one of E_(0,zeta;m),E_(zeta,infinity;m),I_m. The sum of the two endpoint periods is always available, but it is not automatically a decomposition into two single-saddle packets.

## 6. Complete Gaussian packets on regular compact cells

**Theorem CompactCuspThimbleGaussianPackets.** Fix d and a compact subset K of a non-Stokes flow chamber in the upper half-plane. After choosing continuous orientations,

    P_i(m,w)=lambda_i(w)^m*sqrt(2pi/m)*c_i(w)
                                      *[1+O_K(1/m)],
    c_i(w)=1/sqrt(-Phi''(S_i(w),w)) !=0,       (T9)

uniformly for w in K and positive integers m tending to infinity. The square-root branch agrees with the chosen orientation. The same conclusion holds locally on any regular compact set with fixed thimble end labels. No estimate uniform near Stokes walls, the real branch values, w0, or winfinity is asserted.

Proof. In(T4), S=S_i(u) is analytic and the thimble passes along real u from negative to positive. The local integral is lambda_i^m times

    integral_-r^r exp(-m*u^2/2)*S_i'(u)du.

The constant Taylor term is c_i*sqrt(2pi/m); the odd term cancels on the symmetric interval and the remaining local error is O(lambda_i^m*m^(-3/2)). Truncated Gaussian tails are exponentially small.

Every nonlocal arm has f<=f(S_i)-r^2/2. Away from the critical neighbourhoods and the end caps, |Phi'| has a positive lower bound along it. Indeed an arm has constant arg(g_w)=arg(lambda_i); the other critical point has a separated argument, so the arm cannot accumulate there, and its own critical point is excluded by the phase loss. Section4 controls the end caps. Thus |dS/dt| is bounded after the local cut and the complete tail is at most C*|lambda_i|^m*exp(-m*r^2/2)/m. All constants can be chosen uniformly on K: critical points remain simple and separated from the ends and boundaries, the critical arguments stay separated, and the finite-band flow is continuous. This proves(T9) on the complete thimble.

An equivalent exact Laplace representation is available. Parametrize the two arms by g_w(S_i^+(t))=g_w(S_i^-(t))=lambda_i*exp(-t), t>0, with the signs chosen from the local square root. Then

    P_i(m,w)=lambda_i^m*integral_0^infinity exp(-m*t)
                           *[(S_i^+)'(t)-(S_i^-)'(t)]dt.    (T10)

Near t0 the bracket is sqrt(-2/Phi''(S_i))*t^(-1/2)+O(t^(1/2)); Section4 controls infinity. This is another complete-packet proof of the same Gaussian estimate.

## 7. Stokes continuation and the exact rotated endpoint split

The non-Stokes theorem supplies bases on open dense chambers. At a possible Stokes parameter, do not assert that an actual downward arm still avoids the other critical point. Instead continue a neighbouring chamber's end-relative basis. The two endpoint periods and I_m are holomorphic in the upper half-plane: zeta is holomorphic there, the finite path is entire-integrand transport, and a ray through zeta gives an exponentially convergent infinite path. Their exact endpoint relations and(T8) therefore continue to a Stokes parameter. The continued cycles may cease to be single Gaussian packets there. No Stokes-uniform asymptotic estimate is inferred.

Here is the actual ray rotation. Since the integrand is entire and exp(-mS) decays on every large arc inside the closed wedge, rotate the original positive ray to S=zeta*t, t>0. For positive integer m,

    I_m(w)=zeta^(a*m+1)*w^m*[L_m(zeta)+(-1)^m*R_m(zeta)],
    L_m(zeta)=integral_0^1 t^(a*m)*(1-t^d)^m*exp(-m*zeta*t)dt,
    R_m(zeta)=integral_1^infinity t^(a*m)*(t^d-1)^m*exp(-m*zeta*t)dt. (T11)

Both kernels before their complex Laplace factor are positive on their open real intervals. This is an exact identity, not a model-period approximation. It does not by itself make either Laplace transform zero-free or prove two thimble coefficients are nonzero.

## 8. A regular-cell actual noncollision consequence

**Corollary CompactRegularCuspConsecutiveExclusion.** Fix the integer cusp and a compact set K contained in a non-Stokes chamber of the upper half-plane. For all sufficiently large N, the actual pair A_N,A_(N-1) has no common root x with w=x/N^d in K.

Take a slightly larger compact collar inside that chamber; for large N the path

    w_theta=[N/(N-theta)]^d*w, 0<=theta<=1,

stays in the collar. The exact own-order scaling identifies a common root with I_N(w)=I_n(w_1)=0, n=N-1. If only one coefficient in(T8) is nonzero, (T9) excludes a zero. If both are nonzero, write D=log(lambda_R/lambda_L) and B=log|c_R/c_L|. The zero equation and(T9) imply

    F_m(w):=m*Re(D(w))+B(w)=O_K(1/m).         (T12)

The retained sector theorem gives Re(S_R-S_L)>0. Equation(T5) yields exactly

    dD(w_theta)/dtheta=-(S_R-S_L)/(N-theta),
    d[(N-theta)*Re(D(w_theta))]/dtheta
                    =-Re(D(w_theta))-Re(S_R-S_L).          (T13)

At an upper zero, Re(D(w))=O_K(1/N) by(T12); it stays O_K(1/N) on this short path. Also dB(w_theta)/dtheta=O_K(1/N). Compactness makes Re(S_R-S_L)>=c_K>0. Integrating(T13) consequently gives

    F_n(w_1)-F_N(w)<=-c_K+O_K(1/N),

contradicting the O_K(1/N) values forced by two zeros. This proves the corollary. It is an eventual compact-cell statement, not a global all-order theorem or a bound near Stokes/branch/end parameters.


## 9. Both distinguished saddles are active on the entire co-dominant arc

**Theorem CuspEquimodularArcBothActive.** Let H(w)=log|lambda_R/lambda_L|. The distinguished upper-half-plane zero arc H=0 from the positive cusp w_c to its unique negative-axis crossing is contained in one non-Stokes flow chamber. Both coefficients in(T8) are nonzero throughout that arc and its chamber. This asserts activation on the co-dominant arc, not throughout the parameter plane.

Proof. The retained theorem DistinguishedCuspRadialEquimodularGraph supplies a single smooth radial graph, parameterized by arg(w)=theta in(0,pi), with the stated endpoints. Set S_R-S_L=P+iQ; the retained sector theorem gives P>0.

Fix logarithmic lifts without a winding assumption by using the stationary factor itself:

    Phi_i=log(d)+(a+d)*Log(S_i)-Log(S_i-a)-S_i. (T14)

Here Arg(S_i) is in(-pi/d,0) and Arg(S_i-a) is in(-pi,0), so both principal logarithms are globally holomorphic on the corresponding branches. At the positive cusp both saddles tend to s_c>a, giving D=Phi_R-Phi_L tending to0. At the negative crossing S_L tends to a point in(0,a) from below, so Arg(S_L-a) tends to-pi, whereas S_R>a+d tends to the positive axis. Hence Im(D) tends exactly to-pi there. These explicit lifts eliminate a possible hidden multiple of2pi.

Write t=log|w|. Formula(T5) gives

    H_t=-P/d, H_theta=Q/d,
    dt/dtheta=Q/P on H=0,
    d(Im(D))/dtheta=-(P^2+Q^2)/(d*P)<0.        (T15)

Its endpoint limits are0 and-pi. Thus -pi<Im(D)<0 on the whole arc, so it never meets the possible Stokes set. The arc is connected and lies in one component of its complement.

It remains to prove actual activation in that chamber. At the negative real crossing w_0, the positive intervals0 tozeta andzeta toinfinity are precisely the two descending real thimbles, with their unique maxima S_L and S_R. The critical phase difference is-pi, so they are not connected by a descent trajectory. We justify their upper perturbation, including the boundary issue.

Choose analytic local Morse inverses S_i(u,w), real for real u,w near these two real nondegenerate maxima. At u0,

    S_(i,w)=1/Omega'(S_i)<0.

Shrink the local u interval so Re(S_(i,w)(u,w_0))<0 on its entire closure. For w=w_0+i*eta with small eta>0, its image therefore has Im(S_i(u,w))<0 and remains in the lower wedge. These are the initial local thimble arcs. Each of their four finite continuation segments at w_0 reaches a fixed small neighbourhood of its real endpoint or a fixed large infinity cap. Continuous dependence of the regular finite flow segments preserves these entries for upper w perturbations.

The end neighbourhoods need not lie wholly in the wedge. Around0 and the movingzeta, the full small discs are forward-invariant by V=-S/a+O(S^2) and V=-(S-zeta)+O((S-zeta)^2), uniformly in w near w_0. Intersecting them with the separately invariant wedge still forces convergence to the correct end. The infinity cap is uniformly forward-invariant radially by(T6). Thus the perturbed thimbles keep endpoints0-zeta andzeta-infinity. Their tree path from0 toinfinity uses both edges, proving both coefficients are nonzero near the negative anchor.

Local constancy propagates this conclusion through the one non-Stokes chamber containing the entire connected H=0 arc. This proves the theorem without inferring activation merely from the rank of relative homology.


## 10. Classical context and exact remaining scope

Witten, Analytic Continuation of Chern-Simons Theory, arXiv:1001.2933v4, Sections3.1.1-3.1.5, printed24-32, was inspected as primary context for descent cycles, relative homology and Stokes jumps. The one-dimensional Gaussian thimble framework is classical. The present proof verifies the additional logarithmic zeros, inward wedge boundary and noncompact ends directly; those conditions are not supplied merely by citing a boundaryless polynomial-phase example.

The results above determine the complete original-cycle space and its two available saddle contributions on regular chambers, with actual complete Gaussian errors on compact subsets. Both coefficients are proved active on the entire co-dominant arc, while their activation in the rest of the parameter plane has not been classified. No asymptotic bound is uniform at Stokes walls or as w approaches0,w_c,infinity or a real boundary, and no theorem excludes all possible common-root escape to those regions. These are the remaining gaps before a full fixed-cusp Strong Factorial conclusion.

