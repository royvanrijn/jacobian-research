# Cyclotomic filters, adjacent coprimality, and explicit sharp diagonal quadratics

Input: the all-dimension matching ACTUAL-SUPPORT Strong Factorial theorem in
QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md.

Only its diagonal zero-constant specialization is required below. The stronger
k+ell principal is not needed as a separate premise.

For m>=1 fix a primitive m-th root of unity omega, and put

 g_r(X)=X^2+rX,
 F_(m,r)(X_0,...,X_(m-1))=sum_(i=0)^(m-1) omega^i g_r(X_i).

Let mu_j(r)=Lambda_1(g_r^j), and define factorial formal cumulants by

 E_r(T)=sum_(j>=0) mu_j(r)T^j/j!,
 log E_r(T)=sum_(j>=1) kappa_j(r)T^j/j!.

These are formal series over Q[r], not probability-positivity statements.
Define the compressed moment polynomials

 A_(k,m)(r)=Lambda_m(F_(m,r)^(km))/(km)!, k>=0.             (1)

**Theorem 1 (compressed adjacent coprimality).** For every m>=1,
A_(k,m) is a monic integer polynomial of degree km. For every k>=1,
A_(k,m) and A_(k+1,m) are coprime in Q[r].

**Theorem 2 (constructive sharpness).** The polynomial kappa_m(r)
has degree m with leading coefficient (m-1)!, and has no zero root.
Every root r in Qbar gives a phase F_(m,r) with exactly 2m nonzero
diagonal/linear monomials, whose first 2m-1 factorial moments vanish
and whose 2m-th moment is nonzero. The coefficients lie in a number
field of degree at most m*phi(m). Thus the actual-support 2m window
for diagonal quadratics is sharp by an explicit degree-m univariate
construction. This explicit construction concerns start1; the separately
proved homogeneous existence theorem gives sharpness at every start.

**Corollary 3 (doubling cumulant coprimality).** For every m>=1,
kappa_m(r) and kappa_(2m)(r) are coprime in Q[r].

## 1. The formal cyclotomic filter

The factorial functional factors on disjoint variables. Coefficientwise
formal multinomial expansion therefore gives

 sum_(j>=0) Lambda_m(F_(m,r)^j)T^j/j!
      =product_(i=0)^(m-1) E_r(omega^i T)
      =exp(sum_(j>=1) [sum_i omega^(ij)]kappa_j(r)T^j/j!).

The root-of-unity sum is0 unless m divides j, and then is m. Consequently

 product_i E_r(omega^i T)
      =exp(m*sum_(l>=1) kappa_(lm)(r)T^(lm)/(lm)!).        (2)

All moments whose index is not divisible by m vanish identically.
Writing Z=T^m, the right side of(2) is
exp(m*sum_(l>=1) kappa_(lm)(r)Z^l/(lm)!).
This defines all compressed polynomials(1) without any choice of omega.

## 2. Integrality, degree, and adjacent coprimality

Literal expansion gives

 mu_j(r)/j!=sum_(a=0)^j (j+a)!/[a!(j-a)!] *r^(j-a).       (3)

Each coefficient is an integer: it equals binom(j+a,2a)*(2a)!/a!.
Thus each E_r coefficient lies in Z[r]. The product in(2) has
coefficients in Z[omega][r], hence algebraic-integer coefficients.
Every automorphism of Q(omega)/Q permutes its m factors, so these
coefficients are rational. Rational algebraic integers are integers.
This proves A_(k,m) lies in Z[r].

For the leading r-degree, use only the linear term rX_i in each phase.
Its generating series is

 product_i(1-omega^i rT)^(-1)=1/(1-r^m T^m).

Therefore A_(k,m) has leading term r^(km) and no higher degree.

Suppose two adjacent compressed values A_(k,m)(r),
A_(k+1,m)(r) vanished at some algebraic r. Every moment of F_(m,r)
with index between km and (k+2)m-1 would then vanish: the only
multiples of m in that 2m-index interval are km and (k+1)m.
This phase is nonzero, since its diagonal coefficients omega^i are
nonzero. Its actual number s of nonzero monomials is at most 2m (equal to m
at r0). The imported actual-support theorem excludes an initial
subblock of that length. This contradiction proves no common
algebraic root and hence coprimality in Q[r]. This is a universal proof;
no finite resultant computation is a premise.

## 3. Cumulant degree, zeros, and the sharp witness

Each kappa_j is a polynomial of degree at most j, by the formal log
and the degree bounds for mu_j. Its leading r^j coefficient is the
j-th cumulant of the linear phase X. Since

 sum_(j>=0) Lambda_1(X^j)T^j/j!=1/(1-T),
 log(1/(1-T))=sum_(j>=1) T^j/j,

that leading coefficient is (j-1)!, which is nonzero.
Thus kappa_m has degree exactly m and has roots in Qbar.

Choose any such root r. Equation(2) shows every cumulant of F_(m,r)
of orders1,...,2m-1 is zero: nonmultiples of m disappear by the filter,
and the sole remaining order m disappears because kappa_m(r)=0.
Exponentiating, its first2m-1 actual moments are zero.

If r were0, F would have only m diagonal monomials, and its first m
moments would vanish, contradicting the imported actual-support m-window.
Thus r!=0 and all 2m intended diagonal/linear coefficients are nonzero.
The actual-support 2m-window at start1 then forces the 2m-th moment nonzero.
The root r has algebraic degree at most m. Adjoining omega adds degree
at most phi(m), giving the claimed coefficient-field degree bound.
No real/rational witness or minimal field degree is asserted.

Equation(2) also gives

 A_(1,m)=m*kappa_m/m!,
 A_(2,m)=m*kappa_(2m)/(2m)! + (1/2)(m*kappa_m/m!)^2.

Theorem1 applied at k1 therefore proves Corollary3 directly.
Equivalently the nonzero 2m-th moment at every root of kappa_m proves it.

## 4. Concrete dimension-two and dimension-three examples

For m2, kappa_2(r)=r^2+8r+20. Choose r=-4+2i and set

 F_(2,r)=X^2-Y^2+r(X-Y).

It has actual support4, M1=M2=M3=0 and

                      M4=384(12r+59)!=0.

The odd moments vanish by variable exchange; the second moment is
twice kappa_2. This is an explicit witness over Q(i), with nonsingular
diagonal quadratic Hessian.

For m3,

 kappa_3(r)=2(r^3+18r^2+120r+296).
 Equivalently s=r+6 solves s^3+12s+8=0.

Choose ANY root of that cubic and omega with omega^2+omega+1=0.
Then

 F_(3,r)=g_r(X)+omega*g_r(Y)+omega^2*g_r(Z)

has exactly six nonzero monomials, M1=...=M5=0 and M6!=0.
Its coefficients lie in a number field of degree at most6. Literal
cumulant recurrence reduces its sixth moment to

             M6=11520(621r^2+8778r+34264),

whose nonvanishing follows already from the universal upper theorem.
The displayed finite examples are exact illustrations of the universal
proof, not a substitute for it. The construction gives no SFC counterexample.

## 5. Scope, classical relation, and reverse implications

The root-of-unity filter, formal cumulants and product generating series
are derived here coefficientwise. The diagonal actual-support Strong Factorial window is
a genuine imported mathematical premise. The separate stronger k+ell
principal is not required. No positivity or analytic
moment-generating-function convergence is used.

The adjacent coprimality assertion for this compressed one-parameter
family does not imply arbitrary diagonal SFC: this cyclotomic family
does not cover arbitrary diagonal/linear coefficients. The sharpness
witnesses do not alone prove the universal upper bound. The doubling
cumulant coprimality follows one way from compressed adjacency; no
converse is asserted. Full quadratic SFC, SFC(2) and SFC(3) remain open.

No formal/external assurance or worldwide-priority claim is supplied.
