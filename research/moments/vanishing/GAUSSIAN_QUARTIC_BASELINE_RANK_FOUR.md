# A faithful rank-four marked image with a rank-two Gaussian unit

## Exact result

Use one Witt pair Z,W with E[Z^aW^b]=delta_(a,b)a! and an independent
ordinary standard Gaussian Y. Let

    P=W(1+Z)^2Y,  Phi[A]=E[A exp(tP)],
    O=C[t], R=C[t,t^-1], V=Phi[C[Z,W,Y]], M=O V.

Define f=Phi[1] and k=Phi[1-Z]. Then the generic actual marked image has
rank exactly four, with faithful C(t) basis

    1, f, f', k.                                    (1)

All original periods reduce to this basis with coefficients in R.
The EXACT global polynomial-time image is

    M=(R span_C{1,f,f',k}) intersection C[[t]],       (2)

where R span means the R-linear span in the actual formal-series field.
Equivalently M[1/t] is free over R on (1), and M is its origin-regular
saturation. The unit's cyclic module has rank two. The cyclic module of
k has rank three and is a nonsplit primitive extension. A single
original degree-two marker Phi[Z^2] is cyclic of rank four.

The distinction between M and fixed original markers is strict:

    1 belongs to M but not to V.                     (3)

There is a rational horizontal projection pi to the constant line such that

    pi(V)=t^-1 C[t^-1],  pi(M)=R.                    (4)

An explicit finite Laurent-coordinate algorithm below computes the exact
original marker map and its kernel, degree by degree, with no moment-prefix
rank inference. No additional intrinsic characterization of the whole
vector subspace V by asymptotic conditions is claimed.

All formulas have rational coefficients. The same rank, image and kernel
statements hold after any characteristic-zero constant-field extension:
a hypothetical finite rational-function relation uses only finitely many
constants, whose finitely generated characteristic-zero field embeds in C.

## 1. The actual unit, primitive and scalar differential equation

Direct Wick contraction gives, for every a,b,j>=0,

    Phi[Z^aW^bY^j]
      =sum_(m>=0) (m+b)!/m! binom(2m,m+b-a)
                     E[Y^(m+j)] t^m.               (5)

A binomial coefficient is zero when its lower index is outside [0,2m].
Odd ordinary Gaussian moments vanish and E[Y^(2q)]=(2q-1)!!.
This is an exact original-marker formula for every exponent.

In particular

    f=sum_(q>=0) a_q t^(2q),
    a_q=(2q-1)!!binom(4q,2q)
       =(4q)!/[2^q q!(2q)!].

The exact coefficient ratio is

    a_(q+1)/a_q=2(4q+1)(4q+3)/(q+1).

Thus f is divergent at zero. The ratio proves the all-order equation

    16t^3 f''+(48t^2-1)f'+12t f=0.                  (6)

It is the formal 2F0(1/4,3/4;;32t^2), but no hypergeometric theorem
is needed here. Equation (6) has minimal order two: otherwise f'/f
would be rational and regular at zero, since f is a formal unit there.
The corresponding first-order analytic initial-value equation would
make f convergent, a contradiction.

Put g=Phi[1+Z]. Formula (5) gives

    [t^(2q)]g / [t^(2q)]f=(4q+1)/(2q+1),
    k=2f-g=sum_q a_q t^(2q)/(2q+1).

Consequently

    (t k)'=f,    k'=(f-k)/t.                         (7)

The primitive is formal; no convergent integral at zero is asserted.

## 2. An exact Laurent model retaining every original marker

Condition first on Y. The circular root is z=t(1+z)^2Y. Put

    u=1/(1+z), J=1-2t(1+z)Y=2u-1,
    u=(1+sqrt(1-4tY))/2,    s=(u-u^2)/t=Y.

The square root has constant coefficient one. For every integer n define

    F_n=E_Y[u(t,Y)^n/J].                             (8)

Each F_n is a well-defined formal series: u and J are units at t=0,
and every coefficient is polynomial in Y. These are actual conditional
Wick periods. For n<=0 they are original polynomial base markers;
positive n are legitimate rational base amplitudes because 1+z is a
unit at the selected formal root. More explicitly,

    F_n=sum_q (2q-1)!!binom(4q-n,2q)t^(2q),          (9)

using generalized binomial coefficients when the upper index is negative.

The exact spatial phase and normalized volume are

    phi=-(u-u^2)^2/(2t^2),    eta=-du/t.

Ordinary Y Stein integration gives

    2F_(n+3)=3F_(n+2)-F_(n+1)+n t^2F_(n-1),
                      EVERY integer n.             (10)

Indeed phi_u=(-u+3u^2-2u^3)/t^2. On the selected root, u_Y=-t/J,
so (v_u+phi_u v)/J=-(partial_Y-Y)v/t. Its Gaussian expectation is zero
for every Laurent polynomial v(u); all its t coefficients are polynomials
in Y. Applying v=u^n proves (10).

Positive indices reduce by (10) to F0,F1,F2. Negative indices reduce by
solving for F_(n-1) for n<0, leaving F_-1 as the fourth generator.
All reduction coefficients belong to R.

The original W/Y markers admit a particularly short exact amplitude map.
Define Hermite-raising polynomials

    H_(b,j)(Y)=(Y-partial_Y)^b Y^j.

Then the original monomial Z^aW^bY^j maps to the finite Laurent amplitude

    A_(a,b,j)(u,t)
       =(u^-1-1)^a (u^2/t)^b H_(b,j)((u-u^2)/t),     (11)

and its period is the corresponding linear combination of F_n.
To prove this, ordinary Y Stein integration gives

    t Phi[c F(Z)W^bY^j]
      =Phi[F(Z)W^(b-1)Y^(j+1)]
                       -j Phi[F(Z)W^(b-1)Y^(j-1)],
    c=(1+Z)^2.

Iterate in b and use the conditional circular formula at b=0.
Division by c is allowed coefficientwise: c(0)=1, and a fixed t order
uses only finitely many coefficients of the resulting rational Z germ.
This is an actual Wick calculation, not an alternate localized trace.
Equivalently the fixed-z fibre moments have means

    mean(Y)=s, mean(W)=s/(tc),
    cov(WW)=-1/(t^2c^2), cov(WY)=-1/(tc), cov(YY)=0.

The source factor is exp(sigma s/(tc)+xi s-sigma^2/(2t^2c^2)-sigma xi/(tc)).
In particular every term of (11) has form

    constant * t^-p u^n,    p>=0, n<=2p.             (12)

The Z factor only decreases the u exponent.

## 3. Identifying the four actual periods

From (8), u/J=(1+J)/(2J), so

    F0=f, F_-1=g=2f-k, F1=(1+f)/2.

Differentiating the phase and volume in t, or direct coefficient comparison
in (9), gives

    F2=(1+f)/2-2t^2f-4t^3f'.                         (13)

For example f'=-f/t+(F2-2F3+F4)/t^3; use (10) to obtain (13).
The transformation from (F_-1,F0,F1,F2) to (1,f,f',k) has determinant
2t^3, a unit of R. Equations (10)-(13) prove the claimed upper bound
for every original marked period over R, not merely over C(t).

The constant period is genuinely in M:

    Phi[1-t c'(Z)Y]=1,
    f-t Phi[c'(Z)Y]=1.                              (14)

Conditionally this is simply E[J/J]=1. The amplitude in (14) depends on
t, so (14) alone does not put1 in V. The other three basis elements
are literal original periods: f=Phi[1], f'=Phi[P], k=Phi[1-Z].

## 4. Rational primitives prove actual independence

Write (6) as f''=a f'+b f, where

    a=1/(16t^3)-3/t,    b=-3/(4t^2).

If (U f+V f')'=0, independence of f,f' gives

    U'+bV=0, U+V'+aV=0,

hence E(V)=0, where

    E(V)=V''+(1/(16t^3)-3/t)V'
                 +(15/(4t^2)-3/(16t^4))V.           (15)

A rational solution of E(V)=0 has no poles away from zero: at a nonzero
pole, V'' has a uniquely highest pole order. Thus V is a Laurent polynomial.
If its top exponent is N, the top term of E(V) at infinity has coefficient

    (N-3/2)(N-5/2) V_N t^(N-2),

which is nonzero for every integer N. Therefore V=0. It follows that
1,f,f' are independent over C(t): any representation of a nonzero
rational function by U f+V f' can be divided by that function to give1,
then differentiated as above.

If k=U f+V f'+C with rational coefficients, differentiate t k using
(7). Independence of 1,f,f' gives (tC)'=0 and makes tV solve

    E(tV)=-1.                                      (16)

There is no rational solution of (16). Again it must be Laurent polynomial.
The infinity coefficient above forces its top exponent N<=2. Its lowest
exponent M is then at most two; the lowest term at zero is

    (M-3)V_M t^(M-4)/16,

which cannot vanish or match the constant right side. This contradiction
proves independence of all four periods in (1).

Thus the Laurent quartic model is faithful on its generic four-dimensional
cohomology: (10) is an upper spanning reduction, and the four independent
actual images prove the lower bound. The unit module is rank two and
k generates a rank-three nonsplit extension, with quotient represented
by tk. The full image adds its independent constant line.

A low original marker is cyclic for all four dimensions:

    h2=Phi[Z^2]=(1-f)/(2t^2)+f+8t f'+2k.             (17)

In basis (1,f,f',k), stack the coefficient rows of h2,h2',h2'',h2'''.
Their exact determinant is

    3(60t^2+1)/(64t^9)!=0.

This is an algebraic cyclic-rank certificate, not a moment-prefix test.

## 5. The exact global image and original kernel algorithm

The original Gaussian Stein/divergence identity says

    B_f(h)=(h-h(0)f)/t

preserves V and M. Here is a self-contained reason it preserves V. For
each polynomial Q, a centered Hermite expansion supplies a polynomial
vector field X with Gaussian divergence delta(X)=Q-E[Q]. Indeed every
nonconstant Wick-Hermite polynomial is a Gaussian divergence, by its
usual raising identity; the constant coefficient is its Gaussian mean.
Gaussian Stein integration against exp(tP) then gives

    Phi[Q]-E[Q] f=t Phi[X dot grad P].

Thus B_f(Phi[Q]) is again the period of a fixed original polynomial.
For a polynomial-time multiple a(t)Phi[Q], use

    B_f(a h)=a B_f(h)+h(0)(a-a(0))f/t

to obtain preservation of M as well. In particular M is saturated under division by t
whenever its element has zero initial coefficient. Since M contains
(1), while (11)-(13) put every original period in their R span,
M[1/t] is exactly that R span. If h is in this R span and regular at zero,
choose N so t^N h is in the O span of (1). Repeated origin-saturated
division by t puts h in M. This proves the exact global formula (2).
No generic source-rank guess or local formal-solution selection is used.

Equations (10)-(13) and (11) are also a finite exact kernel algorithm.
Expand any original amplitude into monomials; compute (11); reduce its
finitely many u exponents using (10); change frame by (13). This yields
four Laurent-polynomial coordinates in R. The actual period vanishes
if and only if all four coordinates vanish, by independence. Applying
the same procedure to monomials of spatial degree<=D computes the exact
finite-dimensional original marked image and kernel at that degree.
It uses rational algebra and no truncation of the period series.
Formula (5) is the equivalent direct hypergeometric original-marker map.

For example the nonzero original polynomial

    Q=(1+Z)^2(ZW-Y^2)

has zero image. Its Gaussian divergence vector field is
X=(1+Z)^2(W partial_W-Y partial_Y), which is tangent to P, so Stein
gives Phi[Q]=t Phi[X(P)]=0. Adding 1 corrupts this kernel certificate
and gives image f, a nonzero formal unit.

## 6. An exact original-image projection and1notinV

Let pi be the C(t)-linear projection with

    pi(1)=1, pi(f)=pi(f')=pi(k)=0.

It is horizontal because span(f,f',k) is derivative-stable by (6)-(7).
Put ell_n=pi(F_n). Then ell_-1=ell_0=0, ell_1=ell_2=1/2, and (10)
holds for ell_n too. Positive-index induction gives

    ell_n in C[t], deg ell_n<=floor((n-1)/2), n>=1.

Negative-index reduction gives ell_n in t^-1 C[t^-1] for n<0.
Together with (12), every original monomial therefore satisfies

    pi(Phi[A]) in t^-1 C[t^-1].                      (18)

Indeed if n>=1, n<=2p gives deg(t^-p ell_n)<=-1; n=0 has ell0=0,
and negative n already has strictly negative powers. This proves1notinV.

The projection image in (18) is exactly the whole negative-power space.
The fixed original amplitude c'(Z)Y has period j=(f-1)/t by (14), so
pi(j)=-1/t. Since B_f preserves V and pi(f)=0,

    pi(B_f^n j)=-t^(-n-1),    every n>=0.

These are values of fixed original polynomial markers obtained by the
retained spatial Stein operator. They span t^-1 C[t^-1]. Thus (4) holds;
O multiplication gives pi(M)=R. This projection may have poles at zero
although its actual input period is regular there: the poles cancel
between the four faithful formal-series components. It is not an
additional formal-at-zero Gaussian functional.

The complete literal spatial kernel and filtered images also follow
from the all-m [binomial mixed-kernel theorem](GAUSSIAN_BINOMIAL_MIXED_LITERAL_KERNEL.md).
That theorem classifies fixed original markers; the present lattice (2)
classifies polynomial-time amplitudes.

## Verification scope

The proof gives exact original-marker formulas, recurrence reduction,
rational-primitive nonexistence, a cyclic determinant, the global lattice
and the literal projection image. Finite exact controls replay the
identities and compare original Wick coefficients with Laurent reduction.
Those controls do not infer rank or full-image scope from finite jets.
This is internal work, without independent replay, external review or
formal verification claims.

The exact reduction and lattice implementation is [gaussian_quartic_rank_four.py](gaussian_quartic_rank_four.py),
with [bounded exact checker](verify_gaussian_quartic_rank_four.py) and
[verification receipt](gaussian_quartic_rank_four_verification.json).
