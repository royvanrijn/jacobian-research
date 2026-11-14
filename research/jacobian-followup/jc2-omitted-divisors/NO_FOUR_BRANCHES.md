# No affine fibre of either gap-21 normalization has length four

14 November 2026. PF2OMIT12 closes the forbidden Laurent rank<=2 gate
left in PF2OMIT11. It does not classify all Laurent double/triple values
or supply the additional omitted images.

For every coefficient pair in the complete PF2R7FAM1 charts and every
complex target (U,V),

    deg gcd_t(p(t)-U,q(t)-V)<=3.                      (F1)

The bound includes repeated factors and the nonimmersed parameters.
On the immersed domain this is the normalization branch count. At a
nonimmersed point it is a scheme length and must not be substituted for
the number of branches or for an ambient finite-cover length.
The bound is sharp for the combined coefficient family: the L,d=0
triple values in PF2R7COLL2 have length3. No separate sharpness assertion
is made for K.

## 1. A polynomial scaling of the Laurent family

In chart K put y=x^3 and tau=x*t. For x nonzero, scale target coordinates
by x^6 and x^10. Direct substitution in the canonical PF2R7FAM1 formulas
gives polynomials P_y(tau)=x^6*p(tau/x) and Q_y(tau)=x^10*q(tau/x)
in Q[y,tau]. In particular

    P=tau^6+A*tau^4+y*tau^3+C*tau^2+D*tau,
    A=-(15y^2+8y+4)/4,
    C=(225y^4+240y^3+144y^2+48y+16)/48,
    D=-y*(15y^2+4y+4)/12.

The [certificate](four_branch_certificate.json) records the exact Q
and four integer-primitive equations below. The checker reconstructs
both scaled polynomials from the canonical family, rather than treating
the recorded Q or generators as independent assumptions. No division
by a possible immersion or collision factor occurs. Only x!=0, hence
y!=0, will be used.

## 2. Eliminate a hypothetical quartic common divisor

Suppose P-U and Q-V have a common factor of degree at least4. Over C
choose any monic quartic divisor h, allowing repeated roots. Since P-U
is monic of degree6 and has zero tau^5 coefficient, write

    h=tau^4+L*tau^3+M*tau^2+H1*tau+H0,
    (P-U)/h=tau^2-L*tau+g0.

Matching coefficients in degrees4,3,2 forces

    g0=A-M+L^2,
    H1=y+L*M-L*g0,
    H0=C+L*H1-M*g0.

The coefficient of tau gives

    f0=-L*H0+H1*g0-D=0.                              (F2)

The constant target coordinate U is then -H0*g0. Reduce Q modulo the
monic h. Its tau^3,tau^2,tau coefficients must vanish, since its
remainder is the constant V. Let f1,f2,f3 be those three coefficients.
Clearing rational denominators and taking integer primitive parts gives
four generators F0,...,F3 in Q[L,M,y]. This changes each equation only
by a nonzero rational constant. Their supports have17,42,67,85 terms.
All possible quartic factors, including nonreduced ones, are covered;
there is no assumption that their roots or target values are generic.

## 3. Literal power certificate

The saved rational polynomial multipliers C0,...,C3 satisfy

    y^6=C0*F0+C1*F1+C2*F2+C3*F3.                    (F3)

They have386,254,162,147 terms. The
[replay](verify_no_four_branches.py) decodes their polynomial grammar
strictly and multiplies with separate sparse Fraction arithmetic. All
40,579 term products cancel to the one monomial y^6 with coefficient1.
Neither a printed Groebner basis nor modular unit output is a premise.
The discovery scripts and lift log are optional provenance in tmp; the
saved literal identity replays without Singular.

If the four necessary equations vanished at a complex point, (F3)
would force y=0, excluded on K. Hence no quartic common divisor exists.
This proves (F1) on the whole Laurent chart, including its nonimmersed
parameters. The exponent6 is a certificate exponent, not a fibre length
or a claimed minimal exponent.

For L, PF2R7COLL1 exhausts every affine value. On its immersed strata
the maximum number of branches is3. At a nonimmersed parameter its
distinct singular values are one ordinary A2 cusp, one two-branch A5
tangency and eight nodes. The normalization fibre of the cusp has
length2, while each immersed branch contributes length1. Thus no value
there has length greater than2. This proves (F1) on L as well and,
by PF2R7FAM1, for the complete normalized gap-21 coefficient family.

## 4. Omitted-divisor consequence and remaining scope

For either explicit six-by-six image matrix of PF2OMIT11, its corank
is the parametrization fibre length, so its rank is always at least3.
In particular the proposed K rank<=2 exclusion locus is **empty**;
it removes no parameter and must not remain an open elimination task.

On the immersed r7 domain, the known original omitted prime R still
has e2,f1, and every target value has zero, one, two or three distinct
omitted fold points, of total ambient length0,2,4 or6. The matrix rank
computes which case occurs. A three-point value exhausts degree6, so
every additional omitted image must avoid it. This does not make all
other image matching automatic or exclude any remaining global row.

The Laurent triple-value parameter locus, tangency strata, actual B,H
equations, their intersections and the global source realization remain
unresolved. No Keller map is constructed; unrestricted singleton-zero
properness, whole r7 and JC2 remain open. The assurance is a written
cofactor argument, imported complete L classification and exact rational
identity replay, without formalization or external review.
