# An adapted cubic moment excludes the entire degree-thirty class

PF2RULED186. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

The ENTIRE leading degree D=30 is impossible, with unrestricted source
coefficient degrees and unrestricted actual companion degree.
Together with the preceding degree congruence and lower bounds, this gives

    D>=33, b_field=[C(u,w):C(P,Q)]>=7,
    [C(x,t):C(P,Q)]>=14.                                 (1)

These are necessary bounds in the retained complete quadratic class.
The subsequent [square-deformation theorem](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md)
closes every leading degree in this retained marked quadratic ruling-six
model. Global model/import coverage, higher ruling degrees, trace
perfection, normal-tail termination and HC4/JC2 remain open.

The proof uses the first previously unused inverse-prefix equation.
An adapted cubic weight annihilates its negative half-power exactly.
The same weight has only order-seven growth on the finite root cluster,
while its top degree-twenty-seven character has an order-fourteen
escape pole. No bounded companion search or omitted finite contribution
is used.

## 1. The exact residual degree-thirty problem

Retain the original source, full normalizations and volume form:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H^2+R, H=H3, R=d*z+e*w+f.                         (2)

The functions H and R remain source functions; neither is put in
the actual target field or assumed integral over the target.
PF184--185 prove that a D=30 survivor would have

    M=constant*u^2, h=constant*(u-c)^14, c!=0,
    one actual A support c, slope L=7,
    ord_c d=10, ord_c e=3,
    b_field=8, all FOUR escaping sheets poles.           (3)

The entire local equation, with v=u-c and Y=v^7*w, is

    F0(v,Y)-p*v^14=F2*G,
    G=X^4+2J*X^2+q*X+J^2+Delta,
    ord q=10, ord J>=7, ord Delta=14.                   (4)

Its cubic/linear scale X=v^(10/3)*Z has leading quartic

    Z^4+q0*Z, q0 in C*,
    P(v,w)~v^(-2/3)*tau^2*(Z^4+q0*Z),
    dw=v^(-11/3)*dZ, tau in C*.                         (5)

The two other generic roots stay finite, and P(0,w) is a
nonconstant quadratic. No finite or normalized prime is removed.

Every globally reduced actual companion has degree 3k for odd
k>=3. PF185 already proves k=9 in this residual class. The
moment argument below also explains why there is no exception
at a larger index: its last congruence is excluded uniformly.

## 2. The unused inverse-prefix equation

The actual inverse expansion at U=P^(1/6), together with the
whole Jacobian and cubic Kummer symmetry, gives

    Q-sum_(n=0)^k c_n*P^(n/2)
          -c_-1*P^(-1/2)=O(w^-5),
    c_n,c_-1 in C, c_k!=0.                              (6)

In particular the U^-3 coefficient is a complex constant and
the U^-4 coefficient is zero. This is PF185's exact prefix,
not a target change by a half-power.

Write the source cubic as

    H=a3*w^3+a2*w^2+a1*w+a0,
    a3=h/M, a2=K/M.                                    (7)

At the support a2 and a3 have orders seven and fourteen.
Define the NONZERO rational source coefficient

    A=a2/a3=K/h,
    a=-A/2=-K/(2h).                                     (8)

It has pole order seven. No bound on the global degree of K
or of any other coefficient is assumed.

The adapted cubic weight is

    B(w)=w*(w-a)^2/a = w^3/a-2w^2+a*w.                   (9)

Since deg_w R<=2 and deg_w H=3,

    P^(-1/2)=1/a3*w^-3-a2/a3^2*w^-4+O(w^-5).            (10)

Therefore its coefficient residue with this weight is EXACTLY

    [w^-1] B(w)*P^(-1/2)
       =(1/a3)*[-(a2/a3)/a-2]=0.                        (11)

The equality uses a=-(a2/a3)/2. It does not approximate or
discard a negative-half-power contour.

The polynomial Q and every even half-power also have zero
weighted coefficient residue. Multiplication of the remainder
O(w^-5) by the cubic weight gives O(w^-2), so it too has zero
residue. Define, for positive odd n,

    W_n=[w^-1]B(w)*P^(n/2).                              (12)

Taking the coefficient residue of (6) proves the exact ZERO identity

    sum_(1<=n<=k,n odd) c_n*W_n=0.                       (13)

The SAME constants c_n occur as in the earlier companion
relations. This is an identity in C(u), not a numerical
comparison on selected branches.

## 3. The adapted and quartic centers differ only by a bounded term

Let the scaled source cubic factor at its simple finite root as

    S(v,Y)=v^7*H(v,v^-7*Y)
      =C(v)*(Y-r(v))*[(Y-m_H(v))^2+J_H(v)],
    C(0)!=0, m_H(0)=gamma in C*.                        (14)

Its reduction is C(0)*Y*(Y-gamma)^2.

The constant coefficient of H is regular at c because c is off M.
The linear coefficient of H specializes to C(0)*gamma^2!=0,
while its quadratic and cubic coefficients vanish at c.
Thus H has an analytic FINITE simple root rho(v) near c, and

    r(v)=v^7*rho(v)=O(v^7).                              (15)

This is the required stronger estimate at the scaled root
reducing to zero; merely r=O(v) would not suffice.

Comparison of the quadratic/cubic coefficient ratio in (14)
with (7) gives the exact root-sum identity

    A=-v^-7*(r+2m_H),
    a=v^-7*(m_H+r/2).                                   (16)

Let m_G(v) be the midpoint used to depress the separated
quartic at p=0, and put b(v)=v^-7*m_G(v). The source square
first deforms at order t=10. The monic Hensel factorization
and depression therefore give

    m_G-m_H=O(v^10).                                    (17)

No inverse of a vanishing coefficient occurs in that factor
recursion; its reduced factors are coprime with a constant
nonzero resultant. Combining (15)--(17),

    a-b=v^-7*(m_H-m_G+r/2)=O(1).                         (18)

Both centers have the same leading term gamma*v^-7.
All coefficients of the moment center are evaluated at p=0,
because the moments concern P itself. The parameter in the
earlier coframe calculations has not been specialized there.

## 4. The finite cluster now costs only seven

Use the two-cluster contour decomposition proved in PF184.
For small v in a sector, the two finite roots and four
escaping roots are enclosed by separate circles. Both
counts are EVEN, so the branch of sqrt(P) selected at w
infinity is single-valued outside the two discs.

On a FIXED finite circle, P(v,w) converges uniformly to a
nonconstant quadratic with no zero on that circle.
For each fixed positive odd n, P^(n/2) is bounded there.
Since a has pole order seven,

    B(w)=w^3/a-2w^2+a*w=O(v^-7)

on that finite circle. Its contribution to W_n is therefore

    O(v^-7).                                            (19)

This is substantially smaller than the O(v^-21) bound for
the unadapted third centered moment. Equation (11) and
the changed finite bound are both essential. Neither
follows by ignoring the finite roots.

## 5. The complete escape leading coefficient

Use the quartic center b and put

    w=b+epsilon*Z, epsilon=v^(-11/3).

By (18), w-a=epsilon*Z+O(1), and

    w/a=1+O(v^(10/3))

uniformly on a fixed large Z circle. More precisely,
(b-a)/epsilon=O(v^(11/3)), while epsilon/a is
v^(10/3) times a nonzero unit. Consequently

    B(w)=epsilon^2*(Z^2+o(1)).                           (20)

There is no lower-order or extra leading term from the
bounded center error.

The other separated factor is a nonvanishing analytic unit
on this entire escape disc, with square-root limit tau.
Together with (5), the escape contribution is bounded by

    O(v^(-(n+33)/3)),                                   (21)

and its coefficient at that first possible order is

    tau^n*[Z^-1]Z^2*(Z^4+q0*Z)^(n/2).                   (22)

Expansion at Z infinity gives, whenever n is a multiple of three,

    [Z^-1]Z^2*(Z^4+q0*Z)^(n/2)
      =binom(n/2,(2n+3)/3)*q0^((2n+3)/3)!=0.             (23)

For positive odd n the binomial top is a noninteger
half-integer, so this integer coefficient never vanishes.
Uniform convergence on the escape circle justifies taking
this leading integral. Its four roots include the inner
root whose Z limit is zero; it is not excluded from the
contour.

All W_n are rational in u: positive half-power Laurent
coefficients lie in C(u), and a is rational in u.
Thus they have INTEGER Laurent exponents in v.
Combining (19),(21), every positive odd n satisfies

    pole_order W_n<=max(7,floor((n+33)/3)).              (24)

When the nonzero leading order in (23) exceeds seven,
that order determines the EXACT whole pole, since the
finite contribution is strictly later.

## 6. The actual top degree-twenty-seven term cannot cancel

PF185 proves that any residual D=30 actual companion
would have k=9. At n=9, equations (19)--(23) give

    ord_v W_9=-14,
    leading W_9=tau^9*(9/2048)*q0^7!=0.                 (25)

The finite cluster has pole order at most seven.
Every smaller positive odd index n<=7 has, by (24),
pole order at most thirteen. Since c_9!=0, the top
term c_9*W_9 in (13) has a unique nonzero order-fourteen
pole. It cannot cancel. This contradiction excludes
the remaining D=30 class.

For an independent arbitrary-index reading, (23),(24)
exclude EVERY k congruent to three modulo six, k>=3:
its exact top pole (k+33)/3 is at least twelve and
strictly larger than the whole bound for every smaller
positive odd index. PF185's first two positive-moment
relations exclude the other two odd congruences.
Thus no unproved companion-degree cutoff is hidden
in the use of k=9 above.

## 7. New lower bounds and exact remaining scope

PF173 proves that D is divisible by three, with D>=24
and D!=27. PF180 excludes D=24 and the present theorem
excludes D=30. Therefore

    D>=33.

PF166's strict retained-source bound, with m=deg M>=2,
then gives

    b_field>=ceil((D+3m)/6)>=7.

The original quadratic source extension has degree two,
so its original generic map degree is 2*b_field>=14.
This proves (1).

These bounds belong to the complete marked quadratic
ruling-six obstruction model and its retained actual-Keller
hypotheses. They are intermediate implications: PF187--188 subsequently
exclude the entire retained ruling-six class, so D>=33 is no longer a
surviving frontier in that model. Global model/import coverage, higher
rulings, trace integrality, normal-tail termination and HC4/JC2 remain open.

## 8. Exact evidence and replay

The checker verifies the exact adapted negative-character
cancellation, the cubic root-sum center formula, whole
sextic moments with TWO finite simple roots and a nonzero
bounded center mismatch, the order-seven finite bound,
the nonzero all-index cubic congruence coefficient, and
the unique top pole contradiction. The complete written
contour and prefix argument proves the unrestricted result;
finite controls do not replace it.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_adapted_cubic_moment.py

Essential sources:
[PF166](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md),
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md),
[PF180](RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md),
[PF184](RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md),
and [PF185](RULING_SIX_CUBIC_BALANCE_COMPANION_DEGREE_TWENTY_SEVEN.md).
Earlier source/Jacobian and original-normalization hypotheses
are retained through these inputs. Internal derivation and
replay do not constitute formal verification or external
human review.
