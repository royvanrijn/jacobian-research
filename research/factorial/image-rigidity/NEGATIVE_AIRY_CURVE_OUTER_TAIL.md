# A certified cubic loss on the explicit negative-sector Airy tails

This scratch result bounds only the omitted tails of the explicit canonical
Airy contours. It does not prove an original-period nonlocal remainder or a
numerical Strong Factorial threshold. All parameters below are real unless
stated otherwise.

## 1. Exact curve inequalities

For y>0 set

    x=(y-1)*sqrt((y+2)/(3y)), u=x+i*y,
    R=|u|^2=4y^2/3-1+2/(3y),
    D=-Re(u^3/3+u)
     =2(y-1)^2(2y+1)^2/(9y)*sqrt((y+2)/(3y)).

Direct substitution also gives Im(u^3/3+u)=2/3. Its conjugate is the
lower curve. Define

    P=16y^5+64y^4+52y^3-16y^2-y+20.

The following is an exact rational-function identity:

    9 D^2 R-(R-1)^4 = 4(y-1)^4 P/(27y^3).                 (T1)

P is positive: on 0<y<=1, 20-16y^2-y>=3; on y>=1,
52y^3-16y^2-y>=35y^3. All unmentioned terms are positive. Since D>=0,
T1 implies, whenever R>=9/8,

    D/R^(3/2) >= (1-1/R)^2/3 >= 1/243.                  (T2)

There is also a uniform arc-length estimate on exactly this tail set:

    |du| <= 2 |d sqrt(R)|.                                (T3)

Here are complete polynomial certificates for T3. Put

    B=4y^3-3y+2,
    C=4y^4+8y^3+3y^2+2y+1.

Differentiation gives

    (|u'(y)| / |(sqrt(R))'(y)|)^2
       = C B / ((y+2)(4y^3-1)^2).

R is strictly convex, R(3/8)=139/144, and R(1)=1. Therefore R>=9/8
forces either y<3/8 or y>1. On both ranges R is strictly monotone.
Set

    F=4(y+2)(4y^3-1)^2-C B
     =48y^7+96y^6-24y^4-75y^3+3y+6.

In ascending degree order, F(1+z) has coefficients

    [54,594,2079,3429,3096,1584,432,48].

Also (1+z)^7 F(3z/(8(1+z))) has coefficients

    [6,345/8,531/4,114135/512,110697/512,
     60273/512,260283/8192,394833/131072].

Every displayed coefficient is positive. These substitutions cover y>=1
and 0<=y<=3/8, respectively, proving T3 without numerical sampling.

## 2. Physical phase loss and arc length

Let m>=1, 0<q<=1/8, |theta|<=pi/6, and zeta=-q*exp(i*theta).
Scale the two curves by

    t=sqrt(q)*exp(i*theta/2)*u,
    phi(t)=t^3/3-zeta*t.

Their saddle phases have real parts

    c_+=-(2/3)q^(3/2)sin(3theta/2),
    c_-=-c_+.

The real phase on either curve is its c_+ or c_- minus

    L=q^(3/2)cos(3theta/2)D.

For |t|>=rho=3/8, we have R=|t|^2/q>=9/8. Thus T2 and T3 give

    L >= |t|^3/(243 sqrt(2)) > |t|^3/365,
    |dt| <= 2 |d|t||.                                    (T4)

Each curve's outer set has two monotone radial pieces. Consequently, if
T_m denotes the sum of the absolute integrals on both curves restricted
to |t|>=rho, then

    T_m <= 4 (exp(m c_+)+exp(m c_-))
                * integral_rho^infinity exp(-m r^3/365) dr
        <= 4 (exp(m c_+)+exp(m c_-))
                *365/(3m rho^2)*exp(-m rho^3/365).          (T5)

The second inequality follows from
r^3>=rho^3+3rho^2(r-rho). The exact numerical guards are

    365^2>2*243^2,
    rho^3/365=27/186880>1/7000.                            (T6)

Thus this is an explicit outer physical tail bound with positive
absolute loss; the constants are independent of m, q, and theta.

## 3. Uniform comparison with the joint Airy norm

Write xi=m^(2/3)zeta and M(xi)=|Ai(xi)|+|Ai'(xi)|. We use precisely the
two joint-norm bounds proved in EXPLICIT_AIRY_CONTOUR_ENVELOPE.md:

    M(xi)>1/36000                                      if |xi|<=4;
    M(xi)>=(3/(8sqrt(pi)))*|xi|^(-1/4)
                       *(exp(m c_+)+exp(m c_-))        if |xi|>=4.

Their inputs are the actual DLMF Wronskian, connection identity, and
finite asymptotic-error bounds; no separate lower bound for Ai alone is
used. The constants are uniform up to the sector boundaries.

If |xi|<=4, then |m c_+|,|m c_-|<=16/3 and exp(16/3)<225. T5 gives

    T_m/M(xi) < 56064000000*m^(-1)*exp(-m/7000).

If |xi|>=4, T5, sqrt(pi)<2, and |xi|<=m^(2/3)/8 give

    T_m/M(xi) < 19000*m^(-5/6)*exp(-m/7000).

In particular, on the entire stated negative sector,

    T_m < 10^11*m^(-1/3)*exp(-m/14000)*M(xi).             (T7)

The slack factor two in the exponential is deliberately retained for
later polynomial amplitude factors. No such amplitude is silently
included in T7. For zeta=0 choose the two standard decay rays of Ai;
their outer integral is at most
2/(m rho^2)*exp(-m rho^3/3), which satisfies T7 using M(0)>1/36000.

T7 concerns the canonical Airy integrand on complete outgoing curves.
Applying it to an original period requires the separate, proved
matching and amplitude hypotheses on the local inverse domain; original
connecting contours and original nonlocal pieces remain distinct
obligations.

Reproduce the literal polynomial certificates with
`python3 verify_negative_airy_curve_outer_tail.py`.

Research scratch, 20 November 2026. No external human review, formal
verification, historical-priority, or full-period numerical-cutoff claim.

## 4. First and second polynomial weights

Use exactly that note's hypotheses and notation, and let T_(m,j) be
the summed absolute outer-contour integral with the extra factor |t|^j,
where j=0,1,2. The radial arc estimate gives

    T_(m,j) <=4(exp(m c_+)+exp(m c_-))
                  *integral_rho^infinity r^j exp(-m r^3/365) dr.

For 0<=j<=2 the substitution s=r^3 gives the exact elementary bound

    integral_rho^infinity r^j exp(-m r^3/365) dr
      = (1/3)integral_(rho^3)^infinity
                         s^((j-2)/3) exp(-m s/365) ds
      <=365*rho^(j-2)/(3m)*exp(-m rho^3/365).

Consequently all constants in the original comparison are multiplied
by rho^j<=1. Each T_(m,j) therefore satisfies the same uniform bound

    T_(m,j)<10^11*m^(-1/3)*exp(-m/14000)
                          *(|Ai(m^(2/3)zeta)|+|Ai'(m^(2/3)zeta)|).

At zeta=0 the two straight decay rays give the same conclusion directly
by this substitution with cubic coefficient1/3. In particular truncating
the CANONICAL linear amplitude A+B*t has tail bounded by the right side
times |A|+|B|. Arbitrary original amplitudes outside their inverse domain,
original connecting contours, and original nonlocal contributions are
not covered by this statement.

