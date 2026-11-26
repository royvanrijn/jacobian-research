# A gap-free six-window at the quarter-critical constant radius

This proof excludes the constant radius v(k)=3/4 for every naturally integral binary quadratic whose native divided-factorial residue series is nonpolynomial. It imposes no coefficient-grid or positive-gap hypothesis. Polynomial native residues and other constant radii remain separate; this is not unrestricted quadratic SFC.

Let E be a characteristic-zero field with a real-valued nonarchimedean valuation v extending v_2, v(2)=1. Let

    q=alpha X^2/2+beta XY+gamma Y^2/2+dX+eY

have integral natural coefficients in E. Put R_j=Lambda(q^j)/j!, where Lambda(X^i Y^j)=i!j!. The constant k may belong to any real-valued extension. Write O for its valuation ring, enlarging the valued field algebraically when needed.

**Theorem.** If Rbar(T)=sum_j Rbar_j T^j is nonpolynomial and v(k)=3/4, then q+k has no six consecutive zero factorial moments at any positive start.

The proof imports the complete gap-free profile and arbitrary-start parity argument in [the gap-free profile and one-sixth coefficient-field obstruction](BINARY_QUADRATIC_QUARTER_CRITICAL_ONE_SIXTH_COEFFICIENT_FRONTIER.md). The only new step is a four-row universal native coefficient calculation retaining the unknown 2O corrections. No numerical index scan, generic-rank argument or formal-series evaluation at a unit is used.

## 1. All data forced by a hypothetical block

Assume Lambda((q+k)^n),...,Lambda((q+k)^(n+5)) vanish. The imported theorem allows swapping X,Y and dividing the WHOLE phase by an E-unit so that d=1, and gives

    n=8H+2, H>=0,
    C=R_1=alpha+beta+gamma+1+e, v(C-1)=1/6,
    lambda=k^2/2, U=k^4/8, U=C^4 modulo 2O,               (1)

with C of residue one. It also gives the parameter congruences

    alpha=(C+k)(C-1)^2,
    beta=1+lambda+k,
    gamma=C^2(C+k), e=lambda modulo 2O.                  (2)

These are congruences, not exact replacements for individual small moments. That distinction is essential below.

For s=0,...,7 define the ACTUAL small divided moments

    E_s=Lambda((q+k)^s)/s!=sum_(j=0)^s (k^j/j!)*R_(s-j).

Let h in {0,1} be H modulo two. The imported all-H argument proves the following necessary congruences for the original block:

    E_4-2hU in 4O,
    E_5-2hU(C+k) in 4O,
    E_6 in 4O, E_7 in 4O.                               (3)

It also proves E_2,E_3 in 4O, but those two equations will not be needed. Formula (3) retains the TRUE U; it does not set U=C^4 modulo 4O.

## 2. Exact integral parameters, including the constant

Adjoin a root pi of pi^12=2 and extend the valuation; then v(pi)=1/12. Set

    tau=(C-1)/pi^2.

Thus tau is an integral unit. The relation U=C^4 modulo 2O implies

    k=pi^9*C+2chi                                             (4)

for some integral chi. Here is a valuation proof including the possible ramification. Write k=pi^9*(C+x). Since U has residue one, x has positive valuation. If 0<v(x)<1/4, the term x^4 is uniquely least in

    (C+x)^4-C^4=x^4+4Cx^3+6C^2x^2+4C^3x,

and its valuation is less than one. This contradicts U-C^4 in 2O. Hence v(x)>=1/4; dividing pi^9*x by 2 gives integral chi. No completeness, root count or assumption on the coefficient grid is involved.

Because (2) holds modulo 2O, there are A,B,D in O such that the original natural coefficients are EXACTLY

    C=1+pi^2*tau, k=pi^9*C+2chi, lambda=k^2/2,
    alpha=(C+k)(C-1)^2+2A,
    beta=1+lambda+k+2B,
    e=lambda+2D,
    gamma=C-alpha-beta-1-e.                              (5)

The last equation enforces C=the ACTUAL first native moment exactly. It is not an additional assumption. The imported gamma congruence is compatible with (5); no fourth correction variable is discarded by force. A,B,D and chi are arbitrary integral elements subject only to (5).

## 3. Four universal native rows

Put F_4=E_4-2hU, F_5=E_5-2hU(C+k), F_6=E_6, F_7=E_7. Substitution of (5) into the finite native formulas gives the following identities, with every error in pi^16 O:

    F_4 = 2[ h+chi^4+tau^6+pi^3+(1+pi^2*tau)A ],
    F_5 = 2[ (1+h+chi^4)(1+pi^2*tau)+A ],
    F_6 = 2[ 1+tau^6+(1+pi^2*tau)A ],
    F_7 = 2A                              modulo pi^16 O. (6)

In particular all four F_s/2 are integral. The table keeps the original 2O coefficient errors: A occurs explicitly, while B and D cancel to this precision. Individual moments were never asserted invariant under arbitrary corrections in 2O.

Here are complete finite inputs for verifying (6). The calculation takes place in the polynomial ring over Z with integral indeterminates tau,chi,A,B,D, using pi^12=2 and reducing modulo pi^16. Consequently 4=0 and 2pi^4=0, but 2=pi^12 is retained. The parameters may be computed from (5), or from the following short expansions:

    C=1+pi^2*tau,
    k=pi^9+pi^11*tau+pi^12*chi,
    lambda=pi^6+pi^10*tau^2+pi^12*chi^2,
    alpha=pi^4*tau^2+pi^6*tau^3+pi^13*tau^2
          +pi^15*tau^3+pi^12*A,
    beta=1+lambda+k+pi^12*B,
    e=lambda+pi^12*D,
    gamma=C-alpha-beta-1-e                       modulo pi^16. (7)

The true constant exponential coefficients a_j=k^j/j!, j<=7, become

    a0=1,
    a1=pi^9+pi^11*tau+pi^12*chi,
    a2=pi^6+pi^10*tau^2+pi^12*chi^2,
    a3=pi^15,
    a4=1+pi^8*tau^4+pi^12*(1+chi^4),
    a5=a1,
    a6=a2+pi^14*tau^4,
    a7=pi^15                                  modulo pi^16. (8)

Thus U=3a4 at this precision. In particular the pi^12 chi^4 correction to U is included, not lost by replacing U with C^4.

For r=0,...,7 the exact universal normalized native coefficient is

    R_r=sum_(i+j+g+l+m=r)
      [(2i+j+l)!*(2g+j+m)! /
       (2^(i+g)*i!*j!*g!*l!*m!)]
       *alpha^i*beta^j*gamma^g*e^m.                      (9)

Every displayed factorial ratio is an integer, by the usual divided-block counting. Formula (9) contains 792 terms in total through r=7. Multiply (8) with (9), use E_s=sum_j a_jR_(s-j), subtract the two h-dependent quantities defining F_4,F_5, and collect pi-degrees below16. The resulting four complete rows are exactly (6). These formulas are a finite universal integer-polynomial identity certificate; no value of H beyond its parity or unbounded native moment range is being tested. Odd factorial denominators in (8) are units and are reduced modulo4 after removing their powers of two.

## 4. The contradiction

By (3), F_s/2 belongs to 2O. As 2O is contained in pi^4 O, divide (6) by 2 and reduce in the ring O/pi^4 O. This ring has characteristic two. Equation F_7 first gives A=0. Equation F_6 then gives tau^6=1. Since 1+pi^2*tau is a unit, equation F_5 gives chi^4=1+h. The equation F_4 now reduces to

    pi^3=0 modulo pi^4 O.

This is impossible because v(pi^3)=1/4<1/3=v(pi^4). The hypothetical six-zero block cannot exist. This proves the theorem for all positive starts and either parity of H.

The conclusion only removes this critical constant radius in the nonpolynomial-residue sector. It does not assert a six-window for every constant on every coefficient field, does not handle polynomial residues, and does not infer smaller sparse actual-support windows by deleting monomials. For a phase with all six monomials, it does give the actual-support strong factorial window at this radius. All computations used in the proof are exact finite universal identities, not numerical extrapolation.

## 5. The remaining critical radii: a necessary condition

**Corollary.** Under the theorem's natural integrality hypothesis, any
positive six-zero block with nonpolynomial native residue, at any constant
valuation, must have

    v(k)=1-1/L for a power of two L>=8; in particular v(k)>=7/8.

The [previous constant-shift sector theorem](BINARY_QUADRATIC_TWO_POLE_CONSTANT_SHIFTS_STRONG_FACTORIAL_WINDOWS.md),
Section 2, Corollary 3, proves that such a block must have a completely canceled
single-pole residue and a critical constant radius 1-1/L with L>=4.
The theorem above excludes L=4 without any coefficient-grid hypothesis,
leaving L>=8. This is a necessary restriction and does not assert
existence at any remaining radius. It is not a restriction to a value
group with ramification divisible by eight.
