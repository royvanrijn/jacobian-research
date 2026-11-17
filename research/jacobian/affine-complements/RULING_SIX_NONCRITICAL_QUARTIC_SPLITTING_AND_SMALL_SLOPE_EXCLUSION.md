# Noncritical quartic splitting bounds every small pole slope

PF2RULED179. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

At an actual A supporting base value of slope L, retain the balanced
quadratic remainder of PF178. Put b=ord e, t=L+b and let
T(Y)=Y*(A*Y+B) be its scaled leading remainder, while the cubic
square initial is C*Y*(Y-gamma)^2. All A,B,C,gamma are nonzero.

If A*gamma+B!=0, then NECESSARILY

    3t>4(L+1), equivalently 3b>L+4.                       (1)

Consequently every actual A support with L<=5 must satisfy
A*gamma+B=0. In particular the ENTIRE degree-D=24 class, whose
ledger has one support of slope five, admits only that critical
leading alternative. The companion-dependent noncritical ratio
of PF178 is excluded there. The subsequent [quartic pole-budget theorem](RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md)
excludes the entire degree-24 class. Larger-degree critical balanced
cases and full ruling six remain open.

## 1. Exact source hypotheses and scaling

Use the original fields, full normalizations and regular square
form of PF171, with no coefficient or companion degree bound:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H3^2+R2, R2=d(u)z+e(u)w+f(u),
    H3=(h*w+K)*z+B0*w+e0.                               (2)

No source function H3 or R2 is assumed to belong to the target
field or to be integral over the actual target ring.

At an actual A supporting base value c, put v=u-c, Y=v^L*w.
PF178 proves, with a=ord_c d and b=ord_c e,

    a=b+L, 0<=b<=L-2, L>=3,
    M(c)!=0,
    S(v,Y)=v^L*H3(v,v^-L*Y)
       =C*Y*(Y-gamma)^2+v*(higher polynomial in Y),
    v^(2L)*R2(v,v^-L*Y)=v^t*T(v,Y),
    t=L+b<2L, T(0,Y)=Y*(A*Y+B).                         (3)

All the displayed scaled coefficients are in C[[v]], and A,B,C,
gamma are nonzero. The constant term f-d*N/M in R2 is included
in T and has higher order. The entire generic P equation is

    F(v,Y,p)=F0(v,Y)-p*v^(2L)=0,
    F0=S(v,Y)^2+v^t*T(v,Y).                              (4)

The parameter p first occurs at order 2L. This fact is used on
the whole equation, not just its later quartic factor.

## 2. The separated quartic and its first deformation

The simple root zero of S(0,Y) lifts uniquely to r(v) in C[[v]],
r(0)=0. Factor the cubic exactly as

    S=a3(v)*(Y-r(v))*U(v,Y),

where a3 is a unit, a3(0)=C, and U is monic quadratic reducing
to (Y-gamma)^2. The factors Y-r and U are coprime modulo v.

The factors C^2*Y^2 and (Y-gamma)^4 of F modulo v are
coprime. Polynomial Hensel recursion therefore separates F into
a quadratic factor F2 and a MONIC quartic G, with

    F=F2*G, F2(0,Y)=C^2*Y^2, G(0,Y)=(Y-gamma)^4.          (5)

The leading degree-six coefficient is a unit and is assigned to
F2. All coefficients belong to C[p][[v]]: at each recursion step
the only inverse is the fixed constant resultant of the two
reduced factors. Before order 2L they are independent of p.

We compute the quartic correction to S^2 without assuming that
the quadratic factor stays fixed. Introduce an independent
deformation epsilon in S^2+epsilon*T. At epsilon=0 the factors
are a3^2*(Y-r)^2 and U^2. Differentiating their product and
reducing modulo U^2 gives

    partial_epsilon G
      = T/[a3^2*(Y-r)^2] modulo U^2.                    (6)

The right side means its UNIQUE polynomial remainder of degree
at most three, since Y-r is invertible modulo U^2. At v=0,
write X=Y-gamma. That remainder is

    W0(X)=C^-2*[(A+B/gamma)-B*X/gamma^2
                         +B*X^2/gamma^3-B*X^3/gamma^4].  (7)

Indeed it is the Taylor remainder of C^-2*(A+B/Y) through
degree three at Y=gamma. Formula (6) accounts for every
correction to the other factor F2.

After substituting epsilon=v^t, all nonlinear deformation terms
have order at least 2t>t. The term -p*v^(2L) also has order
strictly greater than t. Thus G is U^2 plus v^t*W0 modulo
v^(t+1), after centering the quadratic at its formal midpoint.

Complete the square in the monic quadratic U and then depress
the actual monic quartic. With a new variable X it has the exact form

    G=X^4+2J*X^2+q1*X+J^2+Delta,
    ord_v J>=1 or J=0.                                  (8)

The translations may have p-dependent coefficients beginning at
order 2L; they do not affect any leading term computed here.
The linear coefficient of a depressed quadratic square vanishes
identically. More explicitly, for a square (X^2+J0)^2 perturbed
by epsilon*(w3*X^3+w2*X^2+w1*X+w0), depression gives to first order

    q1=epsilon*(w1-J0*w3),
    Delta=epsilon*(w0-J0*w2).

Since ord J0>=1, (7) therefore proves

    ord_v q1=t,
    initial(v^-t*q1)=-B/(C^2*gamma^2),
    initial(v^-t*Delta)=(A+B/gamma)/C^2.                 (9)

In the NONCRITICAL case A*gamma+B!=0,

    ord_v Delta=t.                                      (10)

All leading coefficients in (8)--(10) are complex constants.
Possible cancellation in J, including J=0, is retained.

## 3. A noncritical quartic always has a simple early root

Put ell=ord_v J when J!=0. Newton expansion of (8) gives at
least one, and in fact at least two, branches with constant
simple leading data and

    mu=ord_v G_X<=3t/4.                                 (11)

Here are ALL cases, including the exceptional repeated initial.

If 2ell<t, center at either of the exact square roots of -J.
Their order is ell/2. The term Delta of order t dominates q1*X,
whose order there is t+ell/2. The displacement scale is
v^((t-ell)/2), smaller than the center because t>2ell.
The initial displacement equation is a nonzero quadratic plus
the nonzero constant initial of Delta. Its two roots are simple.
For all four branches,

    mu=(t+ell)/2<3t/4.                                  (12)

The square roots used for this description may be taken from the
p=0 coefficients; their later p-dependence is irrelevant to the
constant initial. Equivalently apply Newton recursion directly
to the whole equation (4).

If 2ell>t, OR J=0, put X=v^(t/4)*xi. The initial is

    xi^4+delta, delta=initial(v^-t*Delta)!=0.

All four roots are nonzero and simple. Thus

    mu=3t/4.                                            (13)

If 2ell=t, the same scaling gives

    (xi^2+j0)^2+delta, j0,delta nonzero.                  (14)

A NONZERO repeated root would have xi^2=-j0, contradicting
delta!=0. The only possible repeated root is xi=0, which occurs
exactly when delta=-j0^2. In that exceptional case (14) equals

    xi^2*(xi^2+2j0),

and its two nonzero roots are simple. Otherwise all four roots
are nonzero and simple. At any of these nonzero simple roots
one still has mu=3t/4. This includes every leading discriminant
cancellation; the possibly later double-zero cluster is not discarded.

The derivative of the whole F0 along one of these selected roots
has the SAME order mu, since F=F2*G and F2 is a unit near
Y=gamma. A fixed Puiseux substitution v=s^alpha clears the
displayed rational scales. The selected initial equation is
simple with constant derivative. After scaling, (4) has coefficients
in C[p][[s]]; its p term has strictly later order because t<2L.
Constant-unit recursion therefore constructs the ENTIRE branch
with coefficients in C[p]. There is no inversion of a polynomial
in p. In particular the first coefficient depending on p has
affine response, and the same polynomial-coefficient assertion
holds to every later order.

## 4. Exact coframe order forbids small mu

On such a branch, differentiation of the WHOLE equation (4)
with p varying gives

    Y_p=v^(2L)/F0_Y,
    w_p=v^L/F0_Y,
    omega=dQ=-v^L*dv/(M*F0_Y).                           (15)

Its normalized order after v=s^alpha is exactly

    ord_s omega=alpha*(L+1-mu)-1.                       (16)

M(c) is a unit. If mu<L+1, the positive number
alpha*(L+1-mu) is an integer, so (16) is nonnegative.
Thus Q is finite at the branch. Since the entire branch is in
C[p]((s)), its finite residue q0(p) is a polynomial in p.

The source escapes: Y tends to gamma!=0, w has a pole, and
c is off M. In the original t^2=u source cover t is a unit
and the original source x is polar. The retained whole-graph
argument of PF165 applies: this omitted prime has the entire
closed image Q=q0(P), a forbidden affine-line component of the
nonproperness curve. No alternative normalization or selected
subset of primes is used.

If mu=L+1, equation (16) is EXACTLY -1. The coframe has a
nonzero simple-pole coefficient, impossible for dQ in
characteristic zero: the derivative of a rational Laurent series
has zero residue. This rules out equality as well.

Therefore every selected simple branch in Section 3 must satisfy

    mu>L+1.                                             (17)

Together with (11), this proves the uniform necessary inequality
3t>4(L+1), equivalently 3b>L+4, in the noncritical case.

## 5. The small-slope and degree-24 consequences

PF178 supplies b<=L-2. If L<=5, then

    3b<=3L-6<=L+4,

contradicting (1). Hence every such support is CRITICAL:

    L<=5 => A*gamma+B=0.                                 (18)

This also excludes the alternative B/(A*gamma)=-2k/(k-1)
for every odd companion index k>=3 at those supports; that
ratio is never -1.

PF177's COMPLETE degree-24 ledger gives m=2, one A support
with L=5 and intermediate degree [C(u,w):C(P,Q)]=5.
The allowed local orders remain

    (ord d,ord e)=(5,0),(6,1),(7,2),(8,3),

but all now satisfy the critical alternative (18). The deeper
quadratic cluster can have cancellations and nonconstant first
responses. Its complete degree-24 exclusion uses the additional
pole-budget argument of PF180. The larger-slope balanced branch
and full ruling six remain open.

## 6. Exact controls and assurance

The checker verifies the deformation remainder (6),(7), the
depressed coefficients including the corrections to J, all three
Newton cases and the exceptional double-zero initial, exact whole
coframe scalings on full sextic models, and the integer implications.
Its finite controls supplement the arbitrary-coefficient written
proof and do not turn a bounded search into an unrestricted result.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_noncritical_quartic_splitting.py

Essential retained proofs are
[PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md),
[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md),
[PF177](RULING_SIX_DEGREE24_LEADING_CLASSIFICATION_AND_LINEAR_CONTROLS.md),
and [PF178](RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md).
PF177 is needed only for the degree-24 corollary.
No target landing, trace integrality, formal verification, external
human review or complete ruling-six/HC4/JC2 closure is claimed.
