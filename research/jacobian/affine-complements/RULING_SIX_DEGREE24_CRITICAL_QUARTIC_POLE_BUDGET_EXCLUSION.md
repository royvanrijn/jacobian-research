# The entire degree-24 ruling-six class is excluded

PF2RULED180. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

The entire leading degree D=24 is impossible. The complete retained
ledger, actual affine-response pole bound, companion resonance, and
an exhaustive critical quartic Newton calculation exclude all
coefficient choices and every companion degree in that class.

Together with the earlier degree congruence and degree-27 exclusion,
this gives the new necessary bounds

    D>=30, [C(u,w):C(P,Q)]>=6,
    [C(x,t):C(P,Q)]>=12.                                 (1)

The larger-D bounds are necessary conditions. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the balanced
quadratic branch at larger D in the retained marked quadratic model.
HC4/JC2 remain open.

## 1. The complete degree-24 ledger

Use the original source model, full normalizations and volume
form of PF165--PF179. In particular

    P=H3^2+d*z+e*w+f,
    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    d,e nonzero polynomials in u.                        (2)

PF177 proves that D=24 forces

    m=deg M=2,
    b_field=[C(u,w):C(P,Q)]=5,
    exactly ONE actual A supporting base value c,
    its slope L=5.                                     (3)

The companion and all coefficient degrees remain unrestricted.
PF178 gives the balanced local orders

    a=ord_c d=5+b_e,
    b_e=ord_c e in {0,1,2,3}.                            (4)

Write v=u-c, Y=v^5*w and t=5+b_e. The entire equation is

    F=F0(v,Y)-p*v^10,
    F0=S(v,Y)^2+v^t*T(v,Y),
    S(0,Y)=C*Y*(Y-gamma)^2,
    T(0,Y)=Y*(A*Y+B), A,B,C,gamma in C*.                 (5)

By PF179, the noncritical case is already impossible, so

    A*gamma+B=0.                                        (6)

Separate the monic quartic near Y=gamma, as in PF179, and
depress it. Retain every coefficient and branch:

    G=X^4+2J*X^2+q*X+J^2+Delta.                          (7)

All coefficients lie in C(p)[[v]]. Set

    ell=ord_v J>=1, or ell=infinity when J=0,
    ord_v q=t,
    delta=ord_v Delta>=t+1, or delta=infinity.            (8)

PF179's first-variation calculation proves ord q=t even
in the critical case: its leading coefficient is
-B/(C^2*gamma^2), which is nonzero. Equation (6) makes
the order-t Delta coefficient zero. These assertions do not
require a bound on later coefficients or an assumption that
the separated quadratic factor stays fixed.

## 2. Every one of the four sheets must carry an actual pole

For a normalized geometric generic P point above c put
v=s^alpha. Count its sheet multiplicity as alpha. Over all
four escaping quartic points this total sheet mass is FOUR.
The two remaining generic roots are finite and holomorphic;
PF162/PF165 retain all other finite and infinite points.
There is no pole-carrying base value besides c.

At an actual pole, PF164's STRICT affine-response bound gives

    q_Q<3alpha/2,                                       (9)

where q_Q is its positive integral pole order for Q.
If the pole-carrying sheet mass were at most three, summing (9)
would give

    b_field=sum q_Q<9/2, hence b_field<=4,

contradicting (3). Thus ALL FOUR quartic sheets are poles.
This argument also covers ramified points and geometrically
split generic fibres: an index-alpha point contributes alpha
sheets and one pole order q_Q.

At one of these sheets define

    mu=ord_v F0_Y=ord_v G_X.                             (10)

The separated other factor is a unit at the quartic cluster.
PF179's exact whole coframe calculation gives

    ord_s(dQ)=alpha*(6-mu)-1,
    q_Q=alpha*(mu-6).                                   (11)

Combining positivity and (9) therefore gives, on EVERY sheet,

    6<mu<15/2.                                         (12)

Summing (11) with sheet multiplicity alpha and using total
pole degree b_field=5 gives the exact identity

    sum_(four sheets) mu=4*6+5=29.                       (13)

This sum uses geometric sheets, not four different normalized
points when ramification is present. It makes no identification
of a parameter coefficient with a residue-field generator.

## 3. The first two remainder orders are impossible

For a monic depressed quartic as in (7), let xi_i be its four
roots over an algebraic closure of C(p)((v)). Newton sums give

    sum_i xi_i=0, sum_i xi_i^3=-3q,
    sum_i G_X(xi_i)=-8q.                                (14)

Thus at least one sheet has mu<=ord q=t. If t=5 or 6,
this contradicts (12). Consequently

    b_e in {2,3}, t in {7,8}.                            (15)

The argument uses the exact trace of derivative VALUES, not
a comparison between a trace and the sum of their valuations.
The nonzero right side of (14) guarantees that not all four
values can have valuation greater than t.

## 4. Four Newton facts, including every cancellation used below

All valuations here are in v; they may be rational on individual
roots. Coefficients can belong to C(p), so no constant-parameter
or polynomial-residue assumption is needed in this section.

### Square centers

When J is nonzero, put x_+=sqrt(-J), x_-=-x_+.
At these exact square centers,

    G(x_+)=Delta+q*x_+, G(x_-)=Delta-q*x_+.

Let m0=min(delta,t+ell/2). If these two orders differ, BOTH
centers have residual order m0. If they are equal, at least ONE
center has residual order m0: both leading coefficients cannot
cancel because their difference is 2q*x_+!=0.

Suppose

    m0>2ell, m0<2t-ell.                                 (16)

At every center whose residual order is m0, put
X=x_sign+v^((m0-ell)/2)*Z. The square term has initial
4*x_sign^2 times the displacement squared. Its coefficient
has order ell. The displacement is smaller than the center
by m0>2ell, and q times the displacement is higher by the
second inequality in (16). Cubic and quartic displacement
terms are also higher. Hence the initial equation is a
nonzero quadratic plus a nonzero constant. Its two roots
are simple and satisfy

    mu=(m0+ell)/2.                                      (17)

If delta!=t+ell/2, both centers have residual order m0,
so this produces FOUR roots with that mu. If equality allows one
cancellation, it
still produces AT LEAST TWO roots with that mu. Later behavior
of the other two sheets is bounded by (12) and is never omitted.

### A dominant Delta

If

    delta<2ell, delta<4t/3,                              (18)

including ell=infinity, put X=v^(delta/4)*Z. The initial
is Z^4+Delta_lead, with four nonzero simple roots. All have

    mu=3delta/4.                                        (19)

### The square/Delta boundary

If delta=2ell<4t/3, the same scale has initial

    (Z^2+j0)^2+d0, j0,d0 nonzero.

Its only possible repeated root is zero. When d0=-j0^2
it has a double zero and two simple nonzero roots; otherwise
all four roots are nonzero and simple. In EVERY case there
are at least two sheets with

    mu=3delta/4.                                        (20)

This is the entire exceptional leading-discriminant case.

### Dominant cubic/linear balance

If

    ell>2t/3, delta>4t/3,                               (21)

put X=v^(t/3)*Z. The initial is Z*(Z^3+q_lead).
Its three nonzero simple roots have mu=t.
The remaining root is smaller: its order is greater than t/3,
or it is identically zero. There q is the leading derivative,
so its mu is ALSO t. Thus

    all four sheets have mu=t.                          (22)

This last assertion can also be read from the remaining linear
Newton segment: the constant term J^2+Delta has order greater
than 4t/3, while the q coefficient has order t. Possible
cancellation in that constant term only increases the last
root order and does not change its derivative order t.

These facts are formal Newton calculations with simple initial
roots. A fixed Puiseux extension clears their scales, and
unit-derivative recursion over an algebraic closure of C(p)
constructs the indicated roots.
They hold for arbitrary later coefficients.

## 5. Exhaustion for t=7

Here delta>=8, ell is a positive integer or infinity.

If ell<=3, use the square-center fact with
m0=min(delta,7+ell/2). Its hypotheses hold, and at least two
sheets have

    mu<=(7+3ell/2)/2<=23/4<6.

This contradicts (12).

For ell=4, the possibilities are exactly

| Delta order | Newton fact | Consequence |
|---|---|---|
| delta=8 | square/Delta boundary | at least two mu=6 |
| delta=9 | square centers, with possible cancellation at one sign | at least two mu=13/2 |
| delta>=10 | square centers, unequal residual orders | all four mu=13/2 |

The first row contradicts (12). In the second row the other two
sheets have mu<15/2, so the total is STRICTLY below

    2*(13/2)+2*(15/2)=28<29.

The third row gives total 26. Both contradict (13).

For ell>=5, including infinity, the complete list is

| Delta order | Newton fact | Consequence |
|---|---|---|
| delta=8 | dominant Delta | all four mu=6 |
| delta=9 | dominant Delta | all four mu=27/4, total 27 |
| delta>=10 | cubic/linear balance | all four mu=7, total 28 |

Each row contradicts (12) or (13). Thus t=7 is impossible.

## 6. Exhaustion for t=8

Here delta>=9.

If ell<=2, the square-center fact gives at least two sheets with

    mu<=(8+3ell/2)/2<=11/2<6,

contradicting (12). The remaining finite ell and the infinite
case are exhausted by the following table.

| J order | Delta order | Newton fact | Consequence |
|---|---|---|---|
| ell=3 | delta=9 | square centers | all four mu=6 |
| ell=3 | delta>=10 | square centers | all four mu=25/4, total 25 |
| ell=4 | delta=9 | square centers | all four mu=13/2, total 26 |
| ell=4 | delta=10 | square centers, possible one-sign cancellation | at least two mu=7 |
| ell=4 | delta>=11 | square centers | all four mu=7, total 28 |
| ell=5 | delta=9 | dominant Delta | all four mu=27/4, total 27 |
| ell=5 | delta=10 | square/Delta boundary | at least two mu=15/2 |
| ell=5 | delta>=11 | square centers | all four mu=31/4 |
| ell>=6 or infinity | delta=9 | dominant Delta | all four mu=27/4, total 27 |
| ell>=6 or infinity | delta=10 | dominant Delta | all four mu=15/2 |
| ell>=6 or infinity | delta>=11 | cubic/linear balance | all four mu=8 |

Every row directly contradicts (12) or (13), except the
one-sign cancellation row ell=4,delta=10. In that row,
at least two mu=7 and the remaining two are strictly below
15/2. Consequently their total is STRICTLY below

    2*7+2*(15/2)=29,

again contradicting (13). There is no unexamined continuation
of the canceled cluster: the strict bound (12) already excludes
any continuation that could reach the needed total.

This table even covers hypothetical orders delta>=11, so no
assumption about the exact later p-dependent Delta coefficient
or its possible cancellation is needed. Thus t=8 is impossible.

## 7. Consequence and exact scope

Sections 3,5,6 exclude every local order allowed by (4).
Therefore D=24 is impossible with arbitrary coefficients and
an arbitrary-degree actual companion. PF173 already proves
D is divisible by three, D>=24 and D!=27. Hence D>=30.
PF166's retained infinity/source bound then gives

    b_field>=ceil((D+3m)/6)>=6, m>=2.

The original degree is 2*b_field, proving (1).
The entire critical Newton classification needed here is local
at the single compulsory supporting base; none of the four
leading coefficient configurations of PF177 is exempt.

The remaining larger-D balanced quadratic class is not excluded.
No target-field membership, integrality, trace perfection,
normal-tail termination, HC4 or JC2 is inferred.

## 8. Exact controls and assurance

The checker verifies derivative trace (14), sheet-weight and strict
pole-budget identities, all four Newton mechanisms on full quartics,
the exceptional double-zero case, one-sign residual cancellation
controls, the table's exact inequalities and numerical contradictions,
and the strengthened degree bounds. Finite controls supplement the
unrestricted formal proof; they do not prove it by enumeration of
coefficient degrees.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_degree24_critical_exclusion.py

Essential sources:
[PF164](RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md),
[PF166](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md),
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md),
[PF177](RULING_SIX_DEGREE24_LEADING_CLASSIFICATION_AND_LINEAR_CONTROLS.md),
[PF178](RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md),
and [PF179](RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md).
Internal derivation and replay do not constitute formal verification,
external human review or full ruling-six closure.
