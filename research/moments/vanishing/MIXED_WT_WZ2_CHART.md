# The simultaneous WT and WZ^2 chart

14 November 2026. This follows the [cubic-height and linear-WT exclusions](CUBIC_HEIGHT_AIRY_OBSTRUCTION.md).
The full mixed chart is now excluded by the subsequent [rescaled-period
proof](MIXED_PERIOD_AT_INFINITY.md). This note retains the preceding exact
moment input and linear-reduction stage, together with the sharp
eight-moment homogeneous-cubic spine certificate. The earlier standard-basis
run remains an UNKNOWN computation, and the later proof supplies a sufficient 320-moment gate. No global ternary
normal form, modular promotion or new literature-priority claim is made.

## Normalizing the actual remaining family

Consider

    P=W*(a+bZ+dZ^2+cT)+G(Z,T),       deg G<=3, c*d!=0.

The Gaussian-preserving linear null rotation
T -> T+lambda*Z, W -> W-lambda*T-lambda^2*Z/2, Z -> Z, with lambda=-b/c,
removes b. It preserves degree and changes only G besides that coefficient.
Then a nonzero scalar multiple of P and the orthogonal diagonal change
Z -> lambda*Z, W -> lambda^(-1)*W normalize c=d=1 over C. These operations
preserve all pure-vanishing and GMC questions. The new coefficient of W is
an arbitrary parameter alpha (the original a/c). Thus the full remaining
family, without any assumption on alpha, is

    P=W*(alpha+T+Z^2)+G(Z,T),        deg G<=3.       (1)

The eleven coefficients are alpha and the ten coefficients of G. The
normalization does not remove lower-degree terms of G.

## A formal one-variable reduction

For the unnormalized family, let g=t*(a+b*g+d*g^2+c*T). Conditional GL1
gives E(exp(tP))=E_T(exp(t*G(g,T))/(1-bt-2dt*g)). Put

    g=t*(a+cY),
    T=F_t(Y)=Y-(b/c)*t*(a+cY)-(d/c)*t^2*(a+cY)^2.

Then F_t=Y+O(t), and

    dF_t/dY=1-bt-2dt^2*(a+cY)=1-bt-2dt*g.

The Jacobian cancels the conditional Lagrange determinant. Therefore,
coefficientwise in t,

    E(exp(tP))=E_Y exp(Psi(t,Y)),
    Psi=t*G(t*(a+cY),F_t(Y))+(Y^2-F_t(Y)^2)/2.     (2)

Here Y is one standard real Gaussian and Psi is polynomial, of Y-degree
at most six and t-degree at most seven, with Psi(0,Y)=0.

This is a formal Gaussian change, not an analytic contour substitution.
For F_s(Y)=Y+s*(F_t(Y)-Y), differentiating
F_s' exp(-F_s^2/2) Phi(F_s) gives the Y-derivative of
(F_t-Y) exp(-F_s^2/2) Phi(F_s). Gaussian integration by parts makes its
coefficientwise expectation zero. Each t-coefficient is polynomial in Y,
so evaluation at s=0,1 proves the substitution law even when the zero
section moves. Formula (2) follows with its displayed determinant cancellation.

The previous Airy theorem treats real height three; Psi can have height
six. Applying that theorem to (2) without an additional argument would be
invalid. Likewise the one-variable GMC theorem does not apply to a phase
that is nonlinear in t.

## Complete full-parent moment input

Write

    G=r0+u0*T+v0*T^2+s*T^3+r1*Z+u1*Z*T+v1*Z*T^2
          +r2*Z^2+u2*Z^2*T+r3*Z^3.

The [compiler](compile_mixed_wt_parent.py) derives exact integer moment
polynomials through order twelve in all eleven parameters. It does not
sample alpha. For a moment of order m, if k factors come from W*(alpha+T+Z^2)
and the G factors have total Z-degree z0, the number of WZ^2 factors is
forced to be (k-z0)/2. The remaining WT and alpha*W choices are enumerated,
with the exact Wick factors k! and the even T moment. This is a finite
coefficient formula for the full parent, independently checked against
ordinary polynomial expansion.

The first appearances of r0,r1,r2,r3 in moments 1,2,3,4 are respectively

    1,     2alpha,     6(alpha^2+1),
    24alpha*(alpha^2+3).                            (3)

Thus on alpha*(alpha^2+1)*(alpha^2+3)!=0 the first four equations solve
these four coefficients successively. The eight remaining equations from
moments 5,...,12 are exact polynomials in six variables over Q(alpha).
The parameter strata alpha=0, alpha^2+1=0 and alpha^2+3=0 are excluded
from this solver chart and need separate handling in this algorithm. The
subsequent all-order period proof covers them all.

The exact generic elimination was bounded at sixty seconds and 768 MiB
address space. It reached the completed six-variable input, then timed out
in the standard-basis stage. Its authoritative result is **UNKNOWN**.
There is neither a unit certificate nor a surviving all-order polynomial
from that run. The [raw moment input](mixed_wt_parent_moments_12.json),
[solver input](mixed_wt_parent_12.sing),
[run receipt](mixed_wt_parent_12_exact_run.json) and
[log](mixed_wt_parent_12_exact.log) are retained. A separate linear-stage-only
run exported [the completed reduced input](mixed_wt_parent_12_reduced.ideal)
without rerunning the failed standard-basis calculation. Do not restart that
same calculation with larger limits as if the timeout were a mathematical
survivor. An alternate elimination or a new structural constraint is needed.

## A global control for the fifth equation

The retained reduced fifth moment J5 is linear in u2. Write

    J5=K*u2+F,     K=ku*u1+kv*v1+K0,
    ku=120*((a^2+1)^3-4)/(a*(a^2+1)),
    kv=720*(a^4+a^2+2)/(a^2+1),

where a denotes alpha and K0,F are independent of u2; K0 is independent
of u1,v1 as well. The [exact profile](mixed_wt_moment5_profile.json) retains
K explicitly. With x=a^2,

    (-x^2+4x+13)*(x^2+x+2)+(x-6)*((x+1)^3-4)=44.

Consequently the polynomial parameter functions

    ru=a*(a^2-6)*(a^2+1)/5280,
    rv=(-a^4+4a^2+13)*(a^2+1)/31680

satisfy ru*ku+rv*kv=1. On the entire declared open in (3), the change

    K=ku*u1+kv*v1+K0,      L=-rv*u1+ru*v1

has determinant one, with explicit inverse

    u1=ru*(K-K0)-kv*L,     v1=rv*(K-K0)+ku*L.

It introduces no further exceptional a strata. This gives a new route
beyond the timed-out six-variable basis: on K!=0 eliminate u2 from J5,
leaving five variables for moments 6,...,12. On K=0 the independent fifth
equation F=0 must be retained, and u2 is still a variable. These controls alone exclude neither branch; the later period proof
closes the entire chart. The checker reconstructs the first four eliminations
and J5 from the raw integer moments, matches the retained solver input,
and checks this determinant and inverse exactly over Q(a).

## Exact control reductions on all exceptional a faces

The exceptional full-parent faces are not silently discarded. Their first
few moments give the following separate reductions, verified in the exact
quotient fields, with [rational profiles](mixed_wt_exceptional_quadratic_profiles.json)
and an [a=0 profile](mixed_wt_alpha0_first4_profile.json) retained.

* At a=0, moments 1,2,3 eliminate r0,u1,r2 with constant pivots 1,2,6.
  In particular r0=-v0 and u1=-(15s^2+6s*u0+u0^2+2v0^2)/2.
  Put H=u0+5s. The reduced moment four is independent of r3 and has
  u2 coefficient 72H. Reduced moment five has r3 coefficient 360H.
  Thus H!=0 successively eliminates u2,r3, leaving five free coefficients
  for moments 6,...,12. At H=0 retain moment four. Put
  D=r1+10s*v0+5v1+4: reduced moment five has u2 coefficient 360D,
  and moment six has r3 coefficient 2160(D+1). The D=0 and D=-1
  subfaces must be kept if using those further pivots.
* At a^2+1=0, moments 1,2,3,4 eliminate r0,r1,u2,r3 with pivots
  1,2a,12a,48a, all units in Q[a]/(a^2+1). The fifth reduced moment
  is linear in r2 with coefficient

      K1=240*(6a*s*v0-a*u0*v0+3a*v1+2a+6s^2-u1-2v0^2).

  Its u1 slope is the nonzero constant -240, so K1 can replace u1
  without additional parameter exceptions. Split K1!=0 and K1=0.
* At a^2+3=0, moments 1,2,3 eliminate r0,r1,r2 with pivots 1,2a,-12.
  Put H=u0+2s-a*v0. The fourth reduced moment is independent of r3
  and has u2 coefficient -144H. The fifth has r3 coefficient -720H.
  Thus H!=0 again successively eliminates u2,r3, leaving five free
  coefficients for moments 6,...,12. H=0 remains a separate boundary for this reduction.

These are exact low-moment controls for the entire omitted parameter locus,
not exclusions or all-order solutions by themselves. They retain every
boundary needed by this algorithm; the subsequent period proof closes the
chart without promoting the timeout.

## The homogeneous-cubic spine is completely excluded

**Theorem.** There is no alpha,A,B,C,S in C for which the first eight
positive Gaussian moments of

    P=W*(alpha+T+Z^2)+A*Z^3+B*Z^2*T+C*Z*T^2+S*T^3  (4)

all vanish. Eight is sharp for this family: explicit members have the first
seven moments zero. This closes all alpha, including the exceptional values
of the exact elimination below. It does not set the lower terms of G in
(1) to zero without proof.

For legibility use a for alpha in the remainder of this section.

### The a=0 branch

Direct Wick expansion gives

    M2=15S^2.

Thus S=0 at every complex coefficient point. The next useful equations are

    M4=36*(5C^2+6C+4),
    M5=360*B*(5C+4),
    M6=5400*(2A*(C+1)+B^2).

Put q4=5C^2+6C+4. The identities

    12=5q4-(5C+2)*(5C+4),
    3=q4-(5C+1)*(C+1)

force B=A=0. Then M8=529200*q8 with
q8=33C^4+60C^3+60C^2+40C+16, and

    (6567C^3+7452C^2+3912C+3472)*q4
          +(-995C-514)*q8=5664.

This is a contradiction in characteristic zero. The use of S^2=0 here is
pointwise over C; no assertion about an intermediate nilpotent scheme
length is made.

### The a*(a^2+3)!=0 branch

Moments two, three and four successively force

    C=-2-15S^2/(2a),
    B=3S*(75S^2+16a)/(4a^2),
    A=-3*(780S^4*a^2+3375S^4-84S^2*a^3+900S^2*a
                  +16a^4+48a^2)/(8a^3*(a^2+3)).    (5)

These substitutions use exactly the displayed nonzero factors.

If S=0, (5) gives C=-2, B=0, A=-6/a. The sixth and eighth moments reduce
to nonzero multiples, on the declared open, of

    p6(x)=7x^2+72x-45,
    p8(x)=2x^3-135x^2-160x-45,       x=a^2.

They cannot vanish together, since

    (71198x^2-4754655x-9221545)*p6
       +(-249193x-2742363)*p8=538375860.             (6)

If S!=0 set R=S^2. After (5), divide the odd moments by S, clear only
known nonzero denominators, and remove nonzero rational content. The
resulting F5,F6,F7,F8 are stored in the certificate below. In Q[w,a,R],
with F9=w*a*(a^2+3)*R-1, an exact cofactor identity is

    H5*F5+H6*F6+H7*F7+H8*F8+H9*F9=1.              (7)

The [compressed rational certificate](mixed_wt_spine_prefix8_certificate.json.gz)
contains every coefficient. Its 303 cofactor monomials are independently
multiplied and summed with Python Fraction arithmetic. The checker also
reconstructs F5,...,F8 from the original moments; it does not trust an
unmatched solver basis. Identity (7) excludes every point of this open,
including special nonzero values of a, rather than just the generic point
of Q(a).

### The a^2+3=0 branch

Here a!=0, and moments two and three still give C and B from (5).
After clearing nonzero powers of a, moment four is equivalent to

    S^2*(115S^2+128a)=0.

If S=0, C=-2 and B=0. Moment six gives A=32a/5, after which moment eight
is -601316352/5, a contradiction.

If S!=0, S^2=-128a/115, C=146/23 and B=388aS/23. Moment five becomes

    -1440S*(529A+155545a)/529.

Hence A=-155545a/529. Moment six then equals
-3120931595808a/60835, which is nonzero. These computations take place in
the exact quotient by a^2+3 and S^2+128a/115; both roots and both square-root
choices are covered. Together the three branches exhaust C^5 and prove
the eight-moment exclusion.

## A sharp seven-moment false positive

Choose any complex a with

    7a^4+72a^2-45=0,

and set

    P=W*(a+T+Z^2)-(6/a)*Z^3-2Z*T^2.               (8)

No such a is zero or has a^2+3=0. Every monomial has odd circular weight,
so all odd pure moments vanish. Direct contraction gives M2=M4=0 and

    M6=-1440*(7a^4+72a^2-45)/a,
    M8=-241920*(2a^6-135a^4-160a^2-45)/a^2.

Thus moments one through seven are zero. Equation (6) proves the eighth
is nonzero at every displayed root. This is a five-term cubic finite-prefix
trap, not a GMC counterexample or evidence for all-order vanishing.
It also proves that the eight-moment obstruction for (4) is optimal.

The [checker](verify_mixed_wt_chart.py) verifies (2) against direct Wick
expansion, the full-parent moment compiler, the complete spine identities,
the independently multiplied cofactor certificate, exceptional branches
and the exact seven-moment trap. The full parent and all its special a strata are now excluded by the
[rescaled-period proof](MIXED_PERIOD_AT_INFINITY.md). Its optimal finite gate and the other omitted cubic couplings remain open.
