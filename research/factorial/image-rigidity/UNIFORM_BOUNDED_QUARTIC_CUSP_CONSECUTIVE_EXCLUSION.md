# Uniform eventual consecutive exclusion in bounded quartic-cusp windows

This research note proves an eventual order cutoff uniform in the integer cusp parameter, for each fixed bounded scaled field window. The cutoff is not evaluated. There is no global root localization and no full Strong Factorial assertion for cusp supports.

## 1. Statement and normalized family

For each integer s>=1 put a=s^2,d=2s+1,R=(d+1)/(d-1),delta=2d/(d-1), and w_c=a^d*R^(d+1). Define the actual moment polynomial

    A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k).

Set c_d=-delta/sqrt(R)=-2d/sqrt(d^2-1), so 2<=|c_d|<3. For a candidate common root x and upper order N, define the dimensionless scaled field

    h_N=(N^(3/4)/c_d)*[x/(N^d*w_c)-1].       (U1)

**Theorem UniformBoundedQuarticCuspExclusion.** For every finite H>=0 there exists N_0(H), independent of s, such that for all integers s>=1,N>=N_0(H), the actual pair A_N,A_(N-1) has no common complex root satisfying |h_N|<=H.

The window is relative to w_c and scaled as N^(-3/4). This uniform assertion does not apply to a window or zero index growing with N. No numerical N_0(H) is claimed.

Put nu=1/d in [0,1/3], with nu=0 used only as a limiting analytic parameter. For nu>0 write

    delta=2/(1-nu), p=1+2nu*v/(1-nu), R=(1+nu)/(1-nu),
    U=a*(log(p)-(p-1)),
    G_nu(v,q)=U+log[(p^d+q*R^(d+1))/delta].

At nu=0 these are

    G_0(v,q)=-v^2/2+log[(exp(2v)+q*exp(2))/2]. (U2)

The original real v domain is D_nu=(-sqrt(a),infinity), and D_0=R. On each compact domain with p>0 the phase is jointly analytic through nu=0. For real q>0 all logarithms are real. Changing S=a+sqrt(a)*v shows that the original I_m(w_c*q) is the v-period exp(m*G_nu)dv times

    sqrt(a)*[a^(a+d)*delta*exp(-a)]^m,

which is nonzero. Also I_m(w)=(a*m)!/m^((a+d)*m+1)*A_m(m^d*w), so common actual roots compare q_n=(N/n)^d*q_N, n=N-1.

## 2. Uniform geometry and real tails

Fix Q_0>1. All constants in Sections2-4 may depend on Q_0 but not on nu,q,m. Let 1<=q<=Q_0. The real phase has a unique global maximum v_0=v_0(nu,q) in (0,1], with v_0=1 only when q=1. Indeed

    G_v=[-v+delta*p^d/(p^d+q*R^(d+1))]/p.    (U3)

For v>0 let z_q=(q*R^(d+1)+1)/delta, h(v)=(p^d-1)/delta, and

    Z_nu(v)=p^d/v-h(v).

The exact identity

    Z_nu'(v)=-p^(d-1)*(v-1)^2/v^2            (U4)

uses delta*(1-1/d)=2. Thus Z decreases strictly on (0,delta), with limits infinity at0 and 1/delta at delta; z_q>=Z(1). Equation G_v=0 has its unique solution v_0<=1. Formula(U3) is positive for v<0 and negative for v>=delta, so this is the unique global maximum. The same conclusions hold at nu=0 by the limiting formulas.

There is a uniform c_0>0 with c_0<=v_0<=1. For example R^(d+1)<=16 for d>=3, z_q<=(16Q_0+1)/2, and h(v_0)<=15/2; stationarity v_0(z_q+h)=p^d>=1 gives

    c_0=1/[(16Q_0+1)/2+8].                  (U5)

The unique maximizer is continuous in the compact parameter set [0,1/3]x[1,Q_0]; uniqueness and (U5) suffice even at the degenerate endpoint q=1.

The whole real tails admit uniform control. For v<=-1/2 in the original domain, p<=1 and (U3) gives G_v>=-v/p>=1/2. For v>=6, delta<=3 and p<=1+v give

    G_v<=-v/(2p)<=-3/7.                     (U6)

On the fixed compact interval [-1/2,6], the phase and its limiting version are continuous, and the maximum is unique in [c_0,1]. Therefore points a fixed distance from this interval of possible maximizers have a uniform negative phase gap. Together with (U6) this proves exponential bounds on both whole tails, uniformly in nu,q. Every fixed polynomial moment of |v-v_0| is controlled. Factors |G-G(v_0)|^j are also allowed: use t^j*exp(-t)<=C_j*exp(-t/2) for t>=0, then (U6). This also handles the logarithmic divergence at a finite left endpoint p=0. No compact-only argument is being used for an unbounded original tail.

## 3. Uniform quartic-to-Gaussian concentration

Let b=1-v_0>=0. On the fixed interval [c_0/2,2], the positive factors

    A(v)=v/[p*(z_q+h(v))], B(v)=p^(d-1)/v^2

have common positive lower and finite upper bounds. Since

    G_v=-A(v)*integral_(v_0)^v B(t)*(t-1)^2 dt,

integration gives, for r=v-v_0 with v in that interval,

    c*(b^2*r^2+r^4) <= G(v_0)-G(v_0+r)
                         <= C*(b^2*r^2+r^4). (U7)

One universal algebraic core is

    K_b(r)=b^2*r^2/2-b*r^3/3+r^4/12,
    K_b(r)>= (b^2*r^2+r^4)/48,

because the difference equals

    (r^2/16)*(r-8b/3)^2+(5/144)*b^2*r^2.

An upper bound K_b(r)<=b^2*r^2+r^4 follows from the elementary quadratic bound on its cross term. These estimates hold for either sign of r; both integration orientations preserve the positive bounds on A,B.

Moreover

    q-1<=C*b^3.                             (U8)

Indeed z_q-Z(1) is the integral of B(t)*(1-t)^2 from v_0 to1, whose positive factors are bounded on [c_0,1]. Multiplying by delta/R^(d+1) gives (U8).

For the probability density proportional to exp(m*G_nu(v,q))dv on the whole original domain, write E_(m,nu,q). For all sufficiently large real m,

    E[(v-v_0)^2] <= C/(m*b^2+sqrt(m)).        (U9)

Here is a normalization proof. Put ell=(m*b^2+sqrt(m))^(-1/2). On |r|<=ell, the upper bound in (U7) makes the normalized density at least exp(-2C), so its integral is at least c*ell. The upper bound for the unnormalized second moment is C*ell^3: after r=ell*t, either m*b^2*ell^2>=1/2 gives a Gaussian upper envelope, or m*ell^4>=1/4 gives a quartic upper envelope. Outside [c_0/2,2], Section2 supplies C*exp(-c*m); this is bounded by a multiple of ell^3 because ell>=1/sqrt(2m), b<=1. Divide by the lower normalization bound to get(U9).

This is a real-density concentration estimate and uses no complex saddle expansion.

## 4. A positive-side complex zero-free tube

Write w=t*(1+i*eta), t=w_c*q, 1<=q<=Q_0, eta real, |eta|<=1/2. Set

    r(v)=t/(t+S^d)=q*R^(d+1)/(q*R^(d+1)+p^d).

One has |r'(v)|<=3 throughout the original domain. For p<=1 use p^(d-1)<=1 and q*R^(d+1)>=1. For p>=1 use |r'|=delta*r(1-r)/p<=3/4. The limiting nu=0 derivative satisfies the same bound.

Relative to the positive real period at t,

    I_m(w)/I_m(t) = E[(1+i*eta*r(v))^m].

Remove the nonzero phase exp(i*m*eta*r(v_0)). The logarithmic remainder obeys

    |m*log(1+i*eta*r)-i*m*eta*r|<=m*eta^2.

Using |exp(i*x)-1|<=|x| for real x and (U9), the resulting normalized integral differs from1 by at most

    C*m*|eta|/sqrt(m*b^2+sqrt(m))
                          +m*eta^2*exp(m*eta^2).             (U10)

Consequently I_m(w)!=0 whenever the right side is less than1. This is a complete original-period assertion. It covers precisely the positive real side q>=1 needed below, with a uniform bound independent of the cusp gap.


## 5. Uniform bounded-field root branches

For q=1+c_d*epsilon^3*h, epsilon=m^(-1/4), put A=sqrt(R) and v=1+A*epsilon*u. Remove the nonzero scalar A*epsilon*exp(m*G_nu(1,q)) from the v-period. The resulting J_(nu,epsilon)(h) obeys

    J_(nu,epsilon)(h)
       =Z(h)+epsilon*gamma_nu*Z^(5)(h)+O_H(epsilon^2),
    Z(h)=integral_R exp(-u^4/12+h*u)du,
    gamma_nu=nu/[5*sqrt(1-nu^2)].             (U11)

The constants, the first two epsilon derivatives, and each fixed h derivative are uniform in nu in [0,1/3], for bounded h. Here gamma_0=0. A complete justification follows rather than a pointwise-in-d asymptotic import.

Near v=1 the compact nu-family is jointly holomorphic on one fixed neighbourhood. The quartic jet is -2/R^2, bounded away from zero, and the raw mixed derivative G_vvq is zero. The exact same linear-scaling Taylor calculation as in the fixed-cusp proof gives the first two terms of(U11). Its quartic domination and epsilon-derivative bounds have common constants because every local analytic derivative is bounded on this compact family. On the complete real cycle, G_q is uniformly bounded for q in a fixed complex disc about1: p^d is nonnegative real and |p^d+q*R^(d+1)|>=Re(q)*R^(d+1). Thus perturbing q by O_H(epsilon^3) changes the phase uniformly by that amount. Section2 supplies uniform nonlocal gaps and integrability of the base-phase powers generated by epsilon differentiation. Their normalized contributions and first two derivatives are flat, a power of epsilon^(-1) times exp(-c/epsilon^4), uniformly in nu. This proves(U11) on the complete original cycle. At nu=0 use the limiting real domain and the same bounds.

The primary Newman theorem and heat/ODE argument recalled in BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md give that every zero of Z is imaginary and simple. In any fixed compact field window there are finitely many. Uniform convergence, Rouche, and the implicit-function argument therefore give branches

    h_j(nu,epsilon)=h_j+O_H(epsilon),

exhausting all actual zeros in that window, with uniformly bounded first two epsilon derivatives. A one-sided C^2 epsilon function is extended by its quadratic Taylor jet for this use of the ordinary real-parameter implicit-function theorem. Continuity in nu and uniform estimates suffice; no derivative with respect to nu is needed.

For each such branch define the relative w root

    q_j(m;d)=1+c_d*m^(-3/4)*h_j(nu,m^(-1/4)).

Uniformly in d and in this finite branch set,

    |q_j(m;d)-1|<=C_H*m^(-3/4),
    |partial_m q_j(m;d)|<=C_H*m^(-7/4).       (U12)

The derivative holds with the cusp d fixed while m varies. It follows directly by differentiating the C^2 epsilon branch, not an unspecified remainder.

## 6. Exact hard-region reduction

Use the retained corollary cor:linear-lower-exponent in GAMMA_PRODUCT_ROOT_MODULUS.tex, with Gamma beta=1. It proves the actual first N moment polynomials pairwise coprime whenever a>=d(2N-3), including equality. Its statement and proof were read, including the shape-shift hypothesis and actual Gamma-moment identity.

Only the hard case a<d(2N-3) remains. At the integer cusp a=(d-1)^2/4 this implies

    d^2-(8N-10)*d+1<0,
    d<8N-10,
    d/(N-1)<8.                              (U13)

Put n=N-1,kappa=N/n. Hence

    kappa^d<exp(8)<3^8.

For the rest of the proof fix Q_0=2*3^8 in Sections2-4. For bounded upper field |h_N|<=H, the upper relative root is

    q_N=1+c_d*N^(-3/4)*h_N,

so for all large N its real part is in [1/2,3/2]. The lower relative root is exactly q_n=kappa^d*q_N. Its ratio of imaginary to real parts has absolute value at most 6H*N^(-3/4). These estimates are uniform in the hard cusp parameter d.

## 7. Large radial displacement: the lower moment is zero-free

Choose K>=max(1,12H) and suppose d>=K*N^(1/4). For t=Re(w_n), q=t/w_c, the Bernoulli bound kappa^d-1>=d/N gives

    q-1 >= d/(2N)-3H*N^(-3/4) >= d/(4N)>0.  (U14)

Also q<=Q_0 by(U13), for all large N. Section3 gives b>=c*(d/N)^(1/3). In(U10), m=n and |eta|<=6H*N^(-3/4), so its first term is at most

    C_H*(d/N^(1/4))^(-1/3)<=C_H*K^(-1/3).   (U15)

Its second term tends to zero uniformly, bounded by

    36H^2*N^(-1/2)*exp(36H^2*N^(-1/2)).

Choose K=K(H) sufficiently large that(U15) is less than1/4, and then choose N sufficiently large that the second term is less than1/4. Formula(U10) proves I_n(w_n)!=0. Thus no common actual root can occur in this regime. All constants depend only on H and the single fixed Q_0; none depend on d.

## 8. Bounded radial displacement: finite roots have the same index

Now suppose d<=K(H)*N^(1/4). The exact lower field is

    h_n=kappa^(d-3/4)*h_N
                      +c_d^(-1)*n^(3/4)*(kappa^d-1).       (U16)

For all sufficiently large N, uniformly in this regime, kappa^d<=2 and kappa^d-1<=4d/N. Therefore

    |h_n|<=2H+2K(H)=H_1.

Moreover c_d is real, so

    |Im(h_n)-Im(h_N)|<=4H*K(H)*N^(-3/4).     (U17)

The full field difference need not tend to zero; (U17) is the needed replacement for the fixed-d argument.

Use Section5 on a fixed disc larger than H_1+2, with boundary avoiding model zeros. Its finitely many model zeros are all imaginary and have distinct imaginary parts. Uniform Rouche exhaustion places every upper and lower actual zero in a small disc about one of them. Formula(U17) and the O(epsilon) branch deviations imply that upper and lower zeros belong to the same index j, for all sufficiently large N, uniformly in d. If the compact window contains no model zeros there are no actual zeros there and the conclusion is immediate.

A common root would consequently satisfy q_j(n;d)=kappa^d*q_j(N;d). But(U12) gives

    |kappa^d*q_j(N;d)-q_j(n;d)|
       >=(kappa^d-1)*|q_j(N;d)|
                             -|q_j(N;d)-q_j(n;d)|
       >=d/(2N)-C_(H_1)*N^(-7/4)>0          (U18)

for all sufficiently large N. In the second line use |q_j(N;d)|>=1/2, the mean-value bound from(U12) on [n,N], and kappa^d-1>=d/N. Since d>=3, the final strict inequality has a uniform order threshold. This contradicts the common-root equation.

## 9. Quantifiers and remaining scope

Choose K(H) first using Section7, then H_1=2H+2K(H), then a single N_0(H) larger than the finitely many uniform thresholds in Sections2-8. Soft cases are already excluded by Section6; large and bounded radial displacements exhaust every remaining cusp gap. This proves Theorem UniformBoundedQuarticCuspExclusion.

The proof is non-effective: compactness and finite model-zero separation constants have not been evaluated. It does not imply an all-order SFC theorem at any support, localization of all actual roots, or a conclusion for scaled windows growing with N. It improves the fixed-cusp bounded-window theorem only by making the eventual order threshold independent of s. The positive-side tube uses the complete original positive real integral and a concentration estimate; it is not an unproved transfer of model Lee-Yang zeros to the Gamma family.

