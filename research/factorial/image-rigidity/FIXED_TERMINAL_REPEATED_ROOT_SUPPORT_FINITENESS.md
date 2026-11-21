# Fixed-terminal finiteness of repeated-root ordinary binomial supports

For integers m>=2, a>=0 and d>=1 put

    A_(m;a,d)(x)=sum_(k=0)^m binom(m,k)(am+1)_(dk)x^(m-k).

Every polynomial here is monic of degree m and has nonzero constant
coefficient. A support means the exponent pair {a,a+d}; no coefficient
phase is included in the support.

**Theorem FixedTerminalRepeatedRootSupportFiniteness.** For each fixed
integer m>=2, the set

    R_m={(a,d) in Z_(>=0) x Z_(>=1):
                      gcd_Q(A_(m;a,d), A'_(m;a,d)) != 1}

is finite. Equivalently, only finitely many ordinary binomial supports
have an actual degree-m moment polynomial with a repeated complex root.
The finite set and any containing box may depend on m. No uniform bound
in m, numerical containing box, or emptiness assertion is part of this
theorem.

The proof uses only two elementary coefficient criteria already present
in the repository, the finite-degree Lee--Yang circle law, and the
calculations below. It does not use an eventual adjacent-order theorem,
a support-uniform high-order squarefreeness theorem, HJN zero motion,
or a large-degree zero-distribution limit.

## 1. Exact sources and the finite-degree model

The external input is Kabluchko, *Lee--Yang zeroes of the Curie--Weiss
ferromagnet, unitary Hermite polynomials, and the backward heat flow*,
arXiv:2203.05533v3, 30 December 2024. The actual primary PDF was read at
<https://arxiv.org/pdf/2203.05533v3>, also retained in the repository as

    research/factorial/image-rigidity/uniform_cusp_sources/
        kabluchko_2203_05533v3.pdf.

Its equation (1.5), printed p.2, defines

    H_m(z;t)=sum_(j=0)^m (-1)^(m-j) binom(m,j)
                         exp[-t*j*(m-j)/2] z^j.          (F1)

Lemma 2.1, printed p.5, places all its zeros on the unit circle for
every finite m and t>0. This is a finite-degree statement; its proof
credits the classical Lee--Yang theorem. Equation (1.4), printed p.2,
gives the associated entire trigonometric polynomial

    T_m(theta,t)=(2i)^(-m) exp(t*m^2/8)
                        exp(-i*m*theta/2)H_m(exp(i*theta);t)
                 =exp[-(t/2) partial_theta^2]
                                      (sin(theta/2))^m. (F2)

Thus T_t=-T_(theta theta)/2. The exponential-operator notation in (F2)
is unambiguous because there are only m+1 Fourier frequencies. The
external assertion needed below is the finite circle law, not simplicity
or interlacing. Simplicity is proved next.

The internal coefficient inputs are the actual statements and proofs:

* `lower_exponent_initial_and_boundary.tex`, thm:sparse: for a>=1,d>=2,
  d^2>=4(a+d)log(8m) implies individual squarefreeness through order m.
* `ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex`, thm:large:
  2^d>=4m^2 implies squarefreeness of
  F_(m,d)(t)=sum binom(m,k)(dk)! t^k.

Both are under `research/factorial/image-rigidity/`. Their full finite
prefix or all-order conclusions are not inputs here.

## 2. Simplicity of every positive-time finite unitary Hermite model

Fix m and t0>0. By the circle law, all zeros of the entire function
T_m(.,t) are real whenever t>0: a zero at theta would give
|exp(i*theta)|=1, hence Im theta=0.

Suppose T_m(.,t0) has a zero theta0 of multiplicity r>=2. Its finite
Fourier representation is jointly entire in theta and t. Repeatedly
using T_t=-T_(theta theta)/2 in its two-variable Taylor series gives,
locally uniformly for complex w as epsilon->0 through positive reals,

    epsilon^(-r) T_m(theta0+epsilon*w,t0-epsilon^2)
       -> c Q_r(w),
    c=partial_theta^r T_m(theta0,t0)/r! != 0,
    Q_r(w)=exp[(1/2)partial_w^2]w^r
          =sum_(l=0)^floor(r/2)
                r!/[2^l*l!*(r-2l)!] w^(r-2l).          (F3)

For clarity, every Taylor term with theta-degree j and time-degree l
contains a theta derivative of order j+2l. Such terms vanish if
j+2l<r. At j+2l=r their signs cancel because the time increment is
-epsilon^2 and the heat equation contributes (-1/2)^l. All remaining
terms vanish after division by epsilon^r, uniformly on compact w sets.

The monic polynomial Q_r has coefficient of w^(r-1) equal to zero and
coefficient of w^(r-2) equal to r(r-1)/2>0. It cannot have all roots
real: if its real roots were alpha_1,...,alpha_r, their sum would be
zero, forcing

    sum_(i<j) alpha_i*alpha_j
       =-sum_i alpha_i^2/2 <= 0,

contrary to that positive coefficient. Therefore Q_r has a nonreal
root. Choose a small closed disk around such a root, disjoint from the
real axis, whose boundary is zero-free for Q_r. Rouche applied to (F3)
gives a nonreal zero of T_m(.,t0-epsilon^2) for all sufficiently small
epsilon>0. Since t0-epsilon^2>0, this contradicts the circle law.

Consequently every zero of H_m(.;t) is simple for every t>0. This
argument imports neither HJN motion nor a Hermite zero-location theorem.

In particular the monic polynomials

    C_(m,c)(z)=sum_(k=0)^m binom(m,k)
                  exp[-k(m-k)/(2cm)] z^(m-k)
              =(-1)^m H_m(-z;1/(cm)),  c>0,             (F4)

are simple. Their coefficients are continuous in c. On each compact
interval 0<c0<=c<=c1<infinity their discriminants consequently have
a positive minimum modulus. This uses only a finite fixed degree m.

## 3. Uniform actual coefficients when a/d^2 stays positive and finite

Assume a>=1,d>=2 and write c=a/d^2. Scale by W=(am)^d and then by
R(c)=exp(1/(2c)). The exact monic polynomial is

    P_(a,d)(z)=A_(m;a,d)(W*R(c)*z)/(W*R(c))^m
       =sum_(k=0)^m binom(m,k) exp[L_k-k/(2c)] z^(m-k),
    L_k=sum_(j=1)^(dk) log(1+j/(am)).                   (F5)

Uniformly for c in a fixed compact positive interval and 0<=k<=m,

    L_k=k^2/(2cm)+O_(m,c0)(1/d).                        (F6)

Indeed, for every u>=0,

    0<=u-log(1+u)<=u^2/2.

The sum of the linear terms is exactly

    dk(dk+1)/(2am)=k^2/(2cm)+k/(2cmd),

and the error is at most

    (dk)(dk+1)(2dk+1)/(12a^2m^2)=O_(m,c0)(1/d).

There are only m+1 coefficients. Equations (F5)--(F6) therefore give
coefficientwise, locally uniform convergence on every compact c
interval,

    P_(a,d) -> C_(m,c) as d->infinity.                  (F7)

The limiting model is exactly (F4), with time 1/(cm); it is not the
large-degree limiting zero law. By the positive discriminant minimum
from Section 2 and continuity of the discriminant in the coefficients,
P_(a,d) and hence A_(m;a,d) are squarefree for all sufficiently large
d, uniformly when c lies in that compact interval. This controls every
root, without selecting labels or allowing a root to escape a degree-m
polynomial.

## 4. The centered Gamma limit when d/sqrt(a) tends to zero

Here m remains fixed, but d may vary and may stay bounded. Let S have
Gamma shape am+1 and rate m, and put

    V=(S-a)/sqrt(a), delta=d/sqrt(a),
    h(v)=((1+v/sqrt(a))^d-1)/delta.

The support of V is [-sqrt(a),infinity). Direct expansion of the
Gamma integral gives the EXACT monic normalization

    Q_(a,d)(z)
       =A_(m;a,d)(m^d*a^d*(-1+delta*z))
                                      /(m^d*a^d*delta)^m
       =E[(z+h(V))^m].                                 (F8)

No one-shape replacement between different orders is being made;
there is only one fixed actual shape am+1.

The centered Gamma moment-generating function is

    E exp(tV)=exp(-t*sqrt(a))
                   *(1-t/(m*sqrt(a)))^(-(am+1)).        (F9)

As a->infinity its logarithm converges to t^2/(2m), uniformly on every
fixed complex t disk. For example expanding -log(1-u) gives an error
O_(m,T)(a^(-1/2)) on |t|<=T once |t|/(m*sqrt(a))<=1/2.
The linear residual is t/(m*sqrt(a)), the quadratic coefficient is
1/(2m)+1/(2am^2), and the absolute sum of terms of degree at least
three is O_(m,T)(a^(-1/2)). Cauchy consequently proves convergence of
every fixed moment of V to the corresponding moment of a centered
Gaussian of variance 1/m. It also gives a uniform bound on

    E exp(|V|)<=E exp(V)+E exp(-V).                     (F10)

It remains to control the nonlinear map h uniformly in d. For all
v>=-sqrt(a), convexity and Taylor's integral remainder give

    0<=h(v)-v<=(delta/2)*v^2*exp(delta*v_+),
    v_+=max(v,0).                                      (F11)

To check the full support, h(0)=0,h'(0)=1 and

    h''(v)=delta*(1-1/d)*(1+v/sqrt(a))^(d-2).

On the negative half of the support this is between zero and delta.
On the positive half, log(1+v/sqrt(a))<=v/sqrt(a) bounds it by
delta*exp(delta*v). Thus (F11) is global; it does not extrapolate a
bounded-v phase expansion to the far Gamma tail.

If delta<=1/(2m), (F10)--(F11) show

    E|h(V)-V|^m
       <=(delta/2)^m E[|V|^(2m) exp(m*delta*V_+)]
       <=C_m*delta^m ->0.                              (F12)

The last estimate follows from
|v|^(2m)exp(|v|/2)<=C'_m exp(|v|), and (F10). In particular h(V)
and V have uniformly bounded L^m norms. For each 1<=j<=m, the
identity for a difference of powers and Holder's inequality imply

    |E h(V)^j-E V^j|
       <=j*||h(V)-V||_m
                       *(||h(V)||_m+||V||_m)^(j-1) ->0. (F13)

This statement compares variables on their actual common Gamma space;
it does not assert a common-space coupling to a Gaussian. Combined
with the moment convergence obtained from (F9), it proves uniform
coefficient convergence as delta->0:

    Q_(a,d)(z) -> P_m(z)=E[(z+G/sqrt(m))^m],             (F14)

where G is a standard real Gaussian. Here a->infinity automatically,
because d>=2 and delta->0; the error bounds above are independent of d.

The polynomial P_m is monic and has only simple roots. A Gaussian
integration by parts, or direct coefficient comparison, gives

    P_m''+m*z*P_m'-m^2*P_m=0.                           (F15)

If P_m and P_m' vanished at any complex point, uniqueness for this
linear ODE with entire coefficients would force P_m identically zero.
That contradicts monicity. Discriminant continuity applied to (F14)
now shows that A_(m;a,d) is squarefree whenever delta is sufficiently
small, with the threshold depending only on the fixed m. This includes
both bounded d with a->infinity and unbounded d with a/d^2->infinity.

## 5. The high-gap limit and the two boundaries

For a>=1,d>=2, the old sparse theorem applies directly to

    M_(m;a,d)(t)=sum binom(m,k)(am+dk)! t^k
                =(am)!*t^m*A_(m;a,d)(1/t).

All constant and leading coefficients are nonzero, so reversal
preserves multiplicities. In particular A_(m;a,d) is squarefree if

    J=d^2/(a+d)>=4log(8m).                              (F16)

If d->infinity and a/d^2->0, then

    J=1/[a/d^2+1/d] ->infinity,

and (F16) eventually applies.

When a=0, A_(m;0,d)(x)=x^m F_(m,d)(1/x). Its leading and constant
coefficients are nonzero. The raw large-exponent coefficient theorem
therefore excludes repeated roots once 2^d>=4m^2. Only finitely many
gaps d remain on this boundary; no all-order constant-plus-power
finite certificate is used.

When d=1, set gamma=am+1>0. Direct coefficient comparison gives

    x*A''+(1-gamma-m-x)*A'+m*A=0,
    A(0)=(gamma)_m != 0.                               (F17)

A repeated root would be nonzero and an ordinary point of this ODE.
Uniqueness with zero initial value and derivative again contradicts
monicity. Thus d=1 has no repeated-root support for any a>=0.

## 6. A global contradiction at each fixed terminal

Suppose R_m were infinite. Section 5 removes all d=1 supports and all
but finitely many a=0 supports. Choose an infinite sequence of distinct
remaining supports with a>=1,d>=2.

If d is bounded along an infinite subsequence, pass to a fixed d.
Distinctness forces a->infinity, so delta->0 and Section 4 excludes
repeated roots.

Otherwise pass to d->infinity. After taking a subsequence,
c=a/d^2 converges in the extended half-line [0,infinity]. If its limit
is zero, Section 5 applies. If its limit is positive and finite,
Section 3 applies on a compact positive interval containing the tail.
If its limit is infinity, delta=1/sqrt(c)->0 and Section 4 applies.
Every alternative is impossible. This proves R_m finite.

Equivalently, for each fixed m there exists an unevaluated support box
outside which every A_(m;a,d) is squarefree. Rational squarefreeness
also remains valid over every characteristic-zero extension field by
a rational Bezout identity for A and A'. The theorem gives no bound on
that box as m varies, and it does not by itself give a common high-order
cutoff over all supports.

## Source boundary of this scratch result

The sole external mathematical input is the finite Lee--Yang circle
law at positive time, in the exact normalization (F1). The heat
equation, multiplicity blowup, compact-c Gamma coefficient limit,
global Gamma/Hermite moment limit and their simplicity consequences
are proved in this note. The internal logical inputs are precisely
the raw sparse and raw large-exponent coefficient theorems listed in
Section 1. No old terminal-six/ten computation, a0 full finite prefix,
eventual adjacent-order result, five-region theorem, or HJN motion
theorem is a premise. This is a scratch mathematical proof pending
separate review; no source guard, registry entry or canonical edit is
made here.
