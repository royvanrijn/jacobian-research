# A sharp noncritical pole budget and the complete degree-30 ledger

PF2RULED181. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

At any NONCRITICAL actual A support of slope L, put b_e=ord e.
Its entire local Q pole degree satisfies the unrestricted bound

    beta_c<=2*b_e-5<=2L-9.                               (1)

The proof retains the exceptional quartic double-zero cluster and uses
the actual affine first response there. It strengthens the former
2L-5 budget by at least four at every noncritical support.

Consequently the ENTIRE degree-D=30 noncritical case is excluded.
The complete D=30 numerical ledger has only the three families below,
and every survivor in them must have CRITICAL leading remainder
A*gamma+B=0. The subsequent [parameter-precision theorem](RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md)
excludes both slope-six families, and [the nonpole-packet theorem](RULING_SIX_NONPOLE_PACKETS_AND_SLOPE7_TOP_DEGREE_EXCLUSION.md)
excludes intermediate degree nine. The subsequent [two-moment theorem](RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md)
leaves only intermediate degree eight, and the
[companion-prefix theorem](RULING_SIX_CUBIC_BALANCE_COMPANION_DEGREE_TWENTY_SEVEN.md)
forces reduced companion ruling degree twenty-seven.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class. Trace integrality and HC4/JC2 remain open.

## 1. Actual source inputs and a complete noncritical quartic

Retain the original fields, full normalizations and regular square
form of PF171--PF179:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H3^2+d*z+e*w+f, d,e in C[u] nonzero.               (2)

All coefficient and companion degrees are unrestricted. No lower
source function is put into the target field or assumed integral
over the actual target.

At an actual A support c, put v=u-c, Y=v^L*w,
b_e=ord_c e, t=L+b_e. PF178/PF179 give

    ord_c d=L+b_e, 0<=b_e<=L-2, L>=3, t<2L,
    F=F0(v,Y)-p*v^(2L), F0 in C[[v]][Y],
    F=F2*G, F2 unit near Y=gamma,
    G=X^4+2J*X^2+q*X+J^2+Delta,                         (3)

where X is a depressed translation of Y and

    ord q=t, ord Delta=t, ord J=ell>=1 or infinity        (4)

in the NONCRITICAL case A*gamma+B!=0. The translation's
p-dependent coefficients begin only at order 2L.

PF179 proves that this case necessarily has

    3t>4(L+1).                                          (5)

It also constructs the following COMPLETE alternatives.

- If 2ell<t, all four sheets are simple constant-leading
  Newton roots and their derivative orders mu are
  (t+ell)/2<3t/4.
- If 2ell>t or J=0, all four have mu=3t/4.
- If 2ell=t, the scaled initial is
  (xi^2+j0)^2+delta0 with nonzero j0,delta0.
  Unless delta0=-j0^2, all four roots are nonzero simple
  roots with mu=3t/4.
- In the exceptional case delta0=-j0^2, the initial is
  xi^2*(xi^2+2j0). Two outer sheets are simple and have
  mu=3t/4. The remaining TWO inner sheets have ord X>t/4.

Here mu=ord_v F0_Y=ord_v G_X along a sheet.
The pole order of Q at an index-alpha point is

    q_Q=alpha*(mu-L-1)                                  (6)

when it is positive. At every actual pole PF164 gives

    mu<3L/2.                                            (7)

The simple constant-leading sheets have their entire expansions
in C[p] after a fixed Puiseux substitution. Therefore PF179's
whole-graph argument excludes mu<L+1 there; equality would give
an impossible simple pole of dQ. In an actual survivor these
simple sheets all have mu>L+1 and carry poles. The two inner
sheets in the exceptional case may be finite or pole-carrying.

## 2. Exact Taylor precision at an exceptional inner actual pole

We prove a stricter bound at EVERY inner pole in the exceptional
case. Put ell=t/2. Its X order is greater than t/4. On this
entire inner neighborhood, equation (3) gives

    ord_v G_XX=ell, ord_v G_XXX>t/4,
    ord_v G_X>ell, ord_v G>ell.                          (8)

Indeed G_XX=12X^2+4J, where J has exact order ell,
while G_XXX=24X. For G itself the exceptional leading
cancellation makes J^2+Delta have order greater than t;
the other terms also have order greater than t. For G_X,
every term has order greater than 3t/4. All derivative
factorials are nonzero in characteristic zero.

Since F2 and all its derivatives have nonnegative order, and
F2 is a unit, the entire product F0-p*v^(2L)=F2*G has

    ord_v F0_YY=ell EXACT,
    ord_v F0_YYY>t/4,
    ord_v F0_(Y repeated j times)>=0 for every j>=4.     (9)

These bounds are about the WHOLE polynomial, not a derivative
of a selected quotient or a frozen factor.

At an actual pole let the first parameter coefficient of Y
occur at order R, so

    Y=Y0+v^R*(lambda+higher terms),
    Y0 has only C coefficients,
    R=2L-mu>L/2>t/4.                                   (10)

The equality follows from whole implicit differentiation;
strictness follows from PF164's (7). A fixed Puiseux substitution
makes every order integral. The p-dependent depressed translation
starts at 2L>R, so the truncated center Y0 still has its
X coordinate of order greater than t/4. Thus (9) holds at Y0,
with its quadratic coefficient nonzero.

In the Taylor polynomial F0(v,Y0+v^R*Z), the quadratic term has
EXACT order ell+2R. Every cubic term has strictly greater order,
because

    ord F0_YYY+3R>t/4+3R>ell+2R.

Every term of degree j>=4 also has strictly greater order,
since jR>=4R>ell+2R. No later coefficient of the branch can
remove this exact quadratic term in the first parameter
initial. A lowest nonconstant initial below order 2L would
be a nonzero C-polynomial annihilating the first nonconstant
coefficient lambda, which is transcendental over C. This is
impossible. At order 2L the initial is H(lambda)=p. Thus if
ell+2R<=2L, the exact nonzero quadratic term either occurs
below the parameter order or makes that initial genuinely
quadratic; the strictly later higher terms cannot repair it.

PF159 says that this entire first response is AFFINE at an
actual pole. Therefore

    ell+2R>2L,
    mu<L+ell/2=L+t/4.                                  (11)

This is the needed strict inner bound. It does not infer an
affine response from a quadratic truncation; it applies the
actual polynomial-mate theorem to a full Taylor calculation.
The original normalization and every branch remain retained.

## 3. The complete noncritical pole budget

If all four sheets have the simple leading data in Section 1,
each has mu<=3t/4 and must be pole-carrying. Thus

    beta_c<=4*(3t/4-L-1)
           =3t-4L-4
           <2t-2L-4,                                   (12)

because t<2L.

In the exceptional case the two outer sheets contribute exactly

    2*(3t/4-L-1).

An inner pole has charge per sheet, by (6),(11), STRICTLY less
than t/4-1. This upper bound is positive: (5), L>=3 imply t>4.
The total inner sheet mass is two, so their total contribution
is STRICTLY less than 2*(t/4-1). The same strict inequality
holds if only one or neither inner sheet is a pole, since a
finite sheet contributes zero and the bound is positive.
Consequently

    beta_c<2*(3t/4-L-1)+2*(t/4-1)
           =2t-2L-4=2b_e-4.                            (13)

Since beta_c is an integer, (12),(13) prove (1).
This counts index-alpha points with their alpha sheets; it
does not count Puiseux conjugates as different normalized primes.

## 4. The strengthened global ledger

Let n_A>=1 count DISTINCT actual supporting base values, let
n_nc count the noncritical ones, and put S=sum L_c.
PF164's critical bound and (1) give

    b_field=[C(u,w):C(P,Q)]
        <=2S-5n_A-4n_nc.                                (14)

PF166 retains D>=m+4S and b_field>=ceil((D+3m)/6).
Combining those bounds with (14) gives

    4S>=15n_A+12n_nc+2m.                                (15)

In particular if ALL actual supports are noncritical, then

    4S>=27n_A+2m, S>=8, D>=34.

PF173's degree congruence then gives D>=36. The same
strict infinity bound forces b_field>=7 and original
generic degree at least fourteen in that entire subclass.
This is not asserted for a mixed or wholly critical class.

## 5. The complete integer ledger at D=30

Write H=deg h and m=deg M>=2. Since D=m+2H=30,
m is even and H=15-m/2. The complete ledger gives

    5+m/2<=b_field<=2S-5n_A<=H-5n_A.

Thus

    m+5n_A<=10.

The integers m>=2, n_A>=1, m even force n_A=1 and
m in {2,4}. Hence S is the single slope L. The lower
and upper bounds give EXACTLY the following numerical families:

| m | H | L | b_field |
|---|---|---|---|
| 2 | 14 | 6 | 6 or 7 |
| 2 | 14 | 7 | 6,7,8 or 9 |
| 4 | 13 | 6 | 7 |

Indeed L<=floor(H/2), while
5+m/2<=2L-5 forces the displayed lower slopes.
No coefficient search or bounded ledger is used for completeness.

There is only one support. If it were noncritical,
(1) would give beta_c<=3 for L=6 or beta_c<=5 for L=7.
But beta_c=b_field>=6. Therefore the ENTIRE noncritical
degree-30 case is impossible, and EVERY D=30 survivor obeys

    A*gamma+B=0.                                        (16)

All coefficient degrees, later local coefficients and actual
companion degrees remain unrestricted. The three critical numerical families above are necessary ledger data;
the separate successor theorems exclude the slope-six rows and the
slope-seven degree-nine subcase. No realization is asserted.

## 6. Pole mass and some necessary critical order refinements

This paragraph records direct implications of the strict bound
for future work; it is not a classification of critical coefficients.

For L=6, a pole has q_Q<2alpha. If pole sheet mass were at
most three, the total degree would be less than six. Thus
every displayed L=6 family has ALL FOUR sheets pole-carrying.

For L=7, q_Q<5alpha/2. Pole mass at most three gives total
degree at most seven. Thus all four sheets are forced when
b_field is eight or nine; when b_field is six or seven the
pole mass can only be THREE or FOUR.

Whenever all four sheets are poles,
mu>L+1 on each, while the exact derivative trace
sum G_X(xi_i)=-8q gives some mu<=t=L+b_e.
Therefore b_e=0,1 is impossible. This leaves b_e>=2 in
those strata.

There are two further integral refinements. At L=6,b_field=7,
the positive integer pole orders q_Q<2alpha and total
sheet mass four force exactly ONE index-four pole with order
seven. All four mu then equal 7+7/4=35/4, so t>=9 and

    L=6,b_field=7 => b_e>=3.                             (17)

At L=7,b_field=9, the only index partitions capable of
degree nine are one index-four pole of order nine, or
indices (3,1) with orders (7,2). Every mu is then at least
ten; hence t>=10 and

    L=7,b_field=9 => b_e>=3.                             (18)

These are complete integer partition consequences, independent
of a bound on the companion. The remaining critical cases include

| L | b_field | necessary b_e range |
|---|---|---|
| 6 | 6 | 2 through 4 |
| 6 | 7 | 3 or 4 |
| 7 | 6 or 7 | 0 through 5 |
| 7 | 8 | 2 through 5 |
| 7 | 9 | 3 through 5 |

For the L=7,b_field=6 or7 row, orders b_e=0,1 force pole
mass exactly three, by the derivative trace. The table lists
necessary strata, not sufficient data or actual examples.

## 7. Exact controls and remaining scope

The checker verifies the full inner Taylor order argument,
including the cubic and all higher terms, finite index partitions,
the noncritical budgets and strengthened global inequalities,
and the complete D=30 integer reduction. Whole quartic/sextic
controls include the exceptional double-zero cluster and show
why its later nonlinear response cannot be treated as a
simple-root graph automatically.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_noncritical_budget_degree30.py

Essential written inputs are
[PF159](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md),
[PF164](RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md),
[PF166](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md),
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md),
[PF178](RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md),
[PF179](RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md),
and the original complete source model retained there.
Finite controls do not replace those arbitrary-degree proofs.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the degree-30 cubic-balance class described by PF184--185
and the larger-degree mixed and critical classes in the retained marked
quadratic model. Target trace perfection and HC4/JC2 remain open. Independent
whole-proof replay, formal verification and external human
review are unclaimed.
