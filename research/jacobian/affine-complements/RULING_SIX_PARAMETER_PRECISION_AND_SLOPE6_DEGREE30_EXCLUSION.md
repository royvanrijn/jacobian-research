# Parameter precision excludes every slope-six degree-30 stratum

PF2RULED182. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

Every degree-D=30 stratum with actual A slope L=6 is impossible.
Thus any remaining D=30 survivor has exactly one support and

    deg M=2, deg h=14, L=7,
    b_field=[C(u,w):C(P,Q)] in {6,7,8,9},
    M=constant*u^2, h=constant*(u-c)^14, c!=0.            (1)

All coefficient degrees and actual companion degrees remain unrestricted.
The slope-seven critical class is not excluded by this proof.
The subsequent [nonpole-packet theorem](RULING_SIX_NONPOLE_PACKETS_AND_SLOPE7_TOP_DEGREE_EXCLUSION.md)
also excludes degree nine, leaving intermediate degrees six, seven and eight.

The argument also proves a uniform parameter-precision fact:
in the separated depressed quartic, BOTH its constant coefficient
and the correction Delta have a nonzero affine-p coefficient at
order v^(2L). Their generic v orders therefore cannot exceed 2L.

## 1. Original setup and the quartic parameter coefficient

Keep the original source and full normalizations of PF171--181:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H3^2+d*z+e*w+f,
    d,e,f in C[u].                                      (2)

No source function is promoted to the actual target field or
assumed integral. At an A supporting base value c put
v=u-c, Y=v^L*w. The entire generic P equation has the form

    F=F0(v,Y)-p*v^(2L), F0 in C[[v]][Y],
    F(0,Y)=C^2*Y^2*(Y-gamma)^4, C,gamma in C*.            (3)

PF179 separates it into a quadratic factor F2 and a monic
quartic G. Its coefficients belong to C[p][[v]], all coefficients
below order 2L are independent of p, and

    F2(0,Y)=C^2*Y^2, G(0,Y)=(Y-gamma)^4.                 (4)

This Hensel factorization uses only inverses of a constant
nonzero resultant.

Compare the coefficient of p at the FIRST order where p enters,
namely v^(2L). Differentiating the factor product, or just
comparing this coefficient, gives

    -1=(delta F2)*(Y-gamma)^4+C^2*Y^2*(delta G).

Thus the quartic p coefficient is the unique cubic remainder

    delta G=-1/(C^2*Y^2) modulo (Y-gamma)^4.              (5)

Writing X=Y-gamma, its constant coefficient is

    -1/(C^2*gamma^2)!=0.                                 (6)

The midpoint shift of the unperturbed monic quartic reduces to
Y=gamma. Depressing the actual quartic adds a p-dependent
translation only at order 2L. Since its reduction is X^4,
the effect of that translation at this first order has ONLY
a cubic term, and cannot change (6).

Write the exact depressed quartic as

    G=X^4+2J*X^2+q*X+D0
     =X^4+2J*X^2+q*X+J^2+Delta,
    ord_v J>=1, or J=0.                                 (7)

Since J has positive order, the first p coefficient of J^2
has order strictly greater than 2L. Consequently BOTH
D0 and Delta have coefficient at v^(2L) of the form

    constant - p/(C^2*gamma^2).                          (8)

Their earlier coefficients are independent of p. In particular,
over the generic coefficient field C(p),

    ord_v D0<=2L, ord_v Delta<=2L.                       (9)

This proves a generic parameter statement, not an assertion at
every specialized complex p value. It retains all later coefficients.

## 2. The complete slope-six degree-30 pole budget

PF181 proves that a D=30 slope-six stratum has

    L=6, b_field in {6,7},
    exactly one actual pole-supporting base value c,
    ord_c e in {2,3,4},
    t=6+ord_c e in {8,9,10}.                            (10)

The case b_field=7 already has ord e>=3, but it is harmless
to retain the extra t=8 possibility below.

The support is CRITICAL by PF181. The first-variation formulas
of PF179 give

    G=X^4+2J*X^2+q*X+J^2+Delta,
    ord q=t, ord J=ell>=1 or infinity,
    t+1<=delta=ord Delta<=12.                            (11)

The upper bound is (9), so delta is finite.

Every actual pole has positive integral Q pole order satisfying
q_Q<2alpha, where v=s^alpha. If pole sheet mass were at most
three, its total pole degree would be below six, contradicting
b_field>=6. Hence ALL FOUR quartic sheets carry poles.

Let mu be the order of the whole F0_Y, equivalently of G_X.
The exact coframe and strict affine-response bound give

    7<mu<9 on EVERY sheet,
    sum_(four sheets) mu=4*7+b_field in {34,35}.           (12)

The sum counts index-alpha points with alpha sheet multiplicity.
There is no assumption that the four sheets are different
normalized points or that the generic P curve is geometrically
connected.

An immediate exclusion test will be useful:

    if at least TWO sheets have mu<=8,
    their total plus the other two is STRICTLY below 34. (13)

This contradicts (12) for both possible b_field values.

## 3. Newton facts used in the finite order table

We use the arbitrary-coefficient formal Newton facts proved in
PF180, over an algebraic closure of C(p)((v)).

At the square centers x_+=sqrt(-J), x_-=-x_+, the residuals are

    Delta+q*x_+, Delta-q*x_+.

Put m0=min(delta,t+ell/2). If delta!=t+ell/2,
both centers have residual order m0. If the two orders are
equal, at least one center retains order m0, even if the other
leading coefficient cancels. When

    m0>2ell, m0<2t-ell,                                  (14)

there are at least two simple sheets with

    mu=(m0+ell)/2,

and all four have that order when delta!=t+ell/2.

If delta<2ell and delta<4t/3, all four sheets have
mu=3delta/4. If delta=2ell<4t/3, at least two
have that mu, including the exceptional double-zero initial.

If ell>2t/3 and delta>4t/3, cubic/linear balance gives
all four mu=t.

A boundary needed below is equally elementary. A monic depressed
quartic with NONZERO linear coefficient cannot have all its roots
multiple. If it did, its multiplicities would be four or (2,2).
The missing cubic coefficient makes the first case X^4 and the
second case (X^2-a^2)^2, both with zero linear coefficient.
Thus such a leading quartic always has at least one SIMPLE root.

All these claims concern the whole quartic and arbitrary later
coefficients, not a bounded coefficient search.

## 4. Exhaustion at t=8

Here delta is one of 9,10,11,12. If ell<=3, the
square-center calculation gives at least two sheets with

    mu<=t/2+3ell/4<=25/4<8.

For the remaining orders, the complete table is:

| J order | Delta order | Compulsory derivative data |
|---|---|---|
| ell=4 | delta=9 | all four mu=13/2 |
| ell=4 | delta=10 | at least two mu=7, including one-sign cancellation |
| ell=4 | delta=11 or12 | all four mu=7 |
| ell=5 | delta=9 | all four mu=27/4 |
| ell=5 | delta=10 | at least two mu=15/2, including a double-zero initial |
| ell=5 | delta=11 or12 | all four mu=31/4 |
| ell>=6 or infinity | delta=9 | all four mu=27/4 |
| ell>=6 or infinity | delta=10 | all four mu=15/2 |
| ell>=6 or infinity | delta=11 or12 | all four mu=8 |

Every row supplies the pair in (13). Thus t=8 is impossible.

## 5. Exhaustion at t=9

Here delta is 10,11 or12. If ell<=4, the square-center
fact gives at least two mu<=15/2, so (13) applies.

For ell=5:

| Delta order | Compulsory derivative data |
|---|---|
| delta=10 | at least two mu=15/2 |
| delta=11 | all four mu=8 |
| delta=12 | all four mu=33/4, total 33 |

The first two rows satisfy (13), and the third contradicts (12).

For ell>=6, delta=10 gives all four mu=15/2,
and delta=11 gives all four mu=33/4. Both are excluded
by (12),(13). At delta=12 use X=v^3*Z.
The initial quartic is

    Z^4+2j0*Z^2+q0*Z+j0^2+d0,

where j0 may be zero, and q0 is NONZERO. By the last fact
in Section 3 it has at least one simple root, hence at least
one sheet with mu=9. This contradicts the STRICT upper
bound in (12), regardless of every repeated initial root.
Thus t=9 is impossible.

## 6. The sole apparent t=10 possibility has quadratic response

Now delta is 11 or12. If ell<=4, (14) gives at least two

    mu<=t/2+3ell/4<=8,

which is excluded by (13). If ell>=6, delta=11
gives all four mu=33/4, total 33; delta=12 gives
at least two mu=9 (the square/Delta boundary for ell=6
or dominant Delta for ell>6). These contradict (12).

It remains ell=5. At delta=11 all four mu=8, again
excluded. At delta=12 all four mu=17/2 and total 34.
This is the only valuation row not yet excluded. It could
only belong to b_field=6; its actual first response rules
it out.

Choose the square center x_sign=sqrt(-J(v,0)), with either
sign, and the depressed translation evaluated at p=0.
This gives a p-independent Laurent/Puiseux center Yc with
X order 5/2. The parameter-dependent corrections to the
translation start at 12; to the square center they start
at 12-5/2=19/2. They occur far later than the response
order used below.

At this center the exact WHOLE polynomial has

    ord F0(Yc)>=12, ord F0_Y(Yc)>=10,
    ord F0_YY(Yc)=5,
    ord F0_YYY(Yc)>=5/2,
    all derivatives of degree at least four have order>=0. (15)

To verify (15), use the p=0 factorization. At X=x_sign,
G=Delta+q*x_sign has order at least twelve,
G_X=q has order ten, G_XX=-8J has order five,
and G_XXX=24x_sign has order 5/2. Multiplication by
the other, regular unit factor preserves the stated orders.

Substitute

    Y=Yc+v^(7/2)*Z.

The entire equation F0-p*v^12, after division by v^12,
has initial

    a2*Z^2+c0-p, a2 in C*, c0 in C.                    (16)

The linear term has order at least 10+7/2>12.
The cubic term has order at least 5/2+21/2=13>12;
every higher term is later as well. Thus (16) is a GENUINE
quadratic first parameter response. Translating by a constant
coefficient of Yc at the same order cannot change its degree.

The branch has mu=17/2>7 and is an actual Q pole by (12).
PF159 requires its full first response to be affine.
Equation (16) contradicts that necessary condition.
This excludes the final row and therefore every slope-six
degree-30 stratum.

## 7. Exact remaining degree-30 scope

The only complete D=30 ledger row left from PF181 is

    m=2, H=14, L=7, b_field in {6,7,8,9}.

The central M-root has multiplicity at least two, so
M=constant*u^2. The A root c is off M and consumes
ord_c h=2L=14=deg h. Thus h=constant*(u-c)^14
with c!=0, proving (1). The leading primitive equation
has a positive control for these coefficients, so no
leading-equation contradiction is asserted here.

The remaining slope-seven critical quartics, their pole
mass and later companion-resonance constraints remain open.
No derivative integrality, target landing, trace perfection,
normal-tail termination, HC4 or JC2 closure is inferred.

## 8. Exact controls and assurance

The checker verifies the parameter deformation remainder, depression
and first p-constant coefficient, all Newton mechanisms and boundary
polynomials used in the tables, full sextic quadratic-response
controls for the last row, and the exact ledger consequence.
The uniform theorem is the written proof, with no coefficient
or companion degree cutoff.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_slope6_degree30_exclusion.py

Essential retained inputs:
[PF159](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md),
[PF164](RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md),
[PF179](RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md),
[PF180](RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md),
and [PF181](RULING_SIX_NONCRITICAL_POLE_BUDGET_AND_DEGREE30_LEDGER.md).
No independent whole-proof replay, formal verification or external
human review is claimed.
