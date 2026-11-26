# Strong Factorial windows for every univariate polynomial of degree at most two

This proof uses [Dense univariate quadratic Strong Factorial](DENSE_UNIVARIATE_QUADRATIC_STRONG_FACTORIAL.md), Theorem 3, and the binary divided-power splitting construction in its Section 1. The referenced proof digest is d2661315b9c3dab2873d423d23a0c1eab73465520be4e608ce59cc1d428eb364. It concerns ONE variable, not arbitrary binary quadratics.

Let K be a characteristic-zero field and Lambda(x^j)=j!.
For a nonzero polynomial f let N(f) be its actual number of nonzero monomials.

**Theorem 1 (all univariate quadratics).** If deg(f)<=2 and f!=0,
then for every n>=1 at least one of
Lambda(f^n),...,Lambda(f^(n+N(f)-1)) is nonzero.

The dense three-term case is precisely Theorem 3 of the stated proof.
Conversely, the universal assertion of Theorem 1 implies the dense three-term assertion by restriction to abc!=0.
The remaining strata have the following exact arguments.

## 1. An all-unit version of the same divided-power split

For any m>=1, not necessarily even, put eta_k(Z)=2^k Z^k/k!.
Suppose Q(Z)=sum_(k=0)^m b_k eta_k(Z) with EVERY b_k a 2-adic unit.
The same finite-degree Hensel argument factors Q, up to a constant unit,
into factors of distinct degrees equal to the set bits of m.
Every factor has divided-power integral coefficients, a unit constant,
and a unit leading divided-power coefficient. No odd-coefficient
refinement is asserted in this all-unit variant.

Here are the precise differences from the proved construction.
For e the highest bit, m=e+S, the mod2 low algebra A_e is closed,
eta_e eta_j=eta_(e+j) for j<e, and
Qbar=Qlo+eta_e Qhi. All coordinates of Qhi through degree S are1,
in particular its constant and leading coefficients are1.
Choose Ubar=Qhi and Gbar=eta_e+Qlo/Ubar in A_e.
Both constants are1 and G has leading divided-power coefficient1.
The same linearization has high component V and low component
D Ubar+Glo V. It is invertible, with V=Ehi and
D=(Elo-Glo V)/Ubar. Corrections preserve degrees<e and<=S,
and converge coefficientwise in the fixed finite Z_2 lattices.
The residual factor U still has EVERY divided-power coefficient a unit,
because Ubar=Qhi. Recurse until S=0. This includes e=1 and m odd.
No even-index/odd-index requirement was used in these steps.

For a degree e=2^r factor, all positive ordinary coefficients are divisible
by2, while the top coefficient has valuation s_2(e)=1 and the constant
is a unit. Its normalized reciprocal is monic Eisenstein of degree e.
The resulting ordinary polynomial factors are separable and pairwise
coprime: their roots have the distinct valuations1/e. Consequently
the reciprocal polynomial of Q is squarefree over Q_2.
There is only one factor at each degree, because binary digits are distinct.

## 2. Constant plus square: a squarefree moment polynomial

For m>=1 define
g_m(A)=Lambda((Ax^2+1)^m)=sum_(i=0)^m binom(m,i)(2i)!A^i.
Put s=s_2(m) and P_m(B)=2^s g_m(B/4).
The exact factorial valuation identity gives
v_2([B^i]P_m)=s_2(m-i) at EVERY index i, including the endpoints.
Its leading coefficient is a unit. Make it monic by dividing by that unit
and take the reciprocal Q_m(Z)=Z^m P_m(1/Z).
Its ordinary coefficient u_k has valuation s_2(k), so
b_k=u_k k!/2^k is a unit for every k.
Section 1 applies. Thus g_m is squarefree over Q_2 and therefore over Q.

Direct differentiation gives the polynomial identity
A g_m'(A)=m(g_m(A)-g_(m-1)(A)).
Indeed differentiating the phase gives m Lambda(x^2(Ax^2+1)^(m-1)),
and multiplying by A replaces Ax^2 by (Ax^2+1)-1.
If g_(m-1) and g_m had a common zero, it could not be A=0,
since both constant values are1. The identity would make that zero
a repeated root of g_m, contradicting squarefreeness.
Rational Bezout identities transfer this coprimeness to every K.
Scaling f=ax^2+c by its nonzero constant handles every a,c!=0.
Hence this two-term stratum has its exact length-two window.

## 3. Square plus linear term: the elementary Bessel recurrence

Suppose f=ax^2+bx with ab!=0. Define the classical Bessel polynomials
y_m(z)=sum_(i=0)^m (m+i)!/(2^i(m-i)!i!) z^i.
Native factorial expansion gives
Lambda(f^m)=m! b^m y_m(2a/b).
Their exact recurrence is
y_0=1, y_1=1+z, y_(m+1)=(2m+1)z y_m+y_(m-1) for m>=1.

For completeness the coefficient identity behind the recurrence is
(m+i)(m+i+1)-(m-i)(m-i+1)=2i(2m+1);
multiplying by (m+i-1)!/[2^i(m+1-i)!i!] verifies every coefficient,
with out-of-range coefficients interpreted as zero.
A simultaneous zero of consecutive y_m,y_(m+1) propagates backward
through the recurrence to y_0=0, impossible.
Thus this stratum also has its exact length-two window, over every K.
No Bessel irreducibility theorem is used.

## 4. Affine and monomial strata

For f=bx+c with bc!=0, set E_m(z)=sum_(j=0)^m z^j/j!.
Then Lambda(f^m)=m!b^m E_m(c/b) and
E_(m+1)-E_m=z^(m+1)/(m+1)!.
Consecutive vanishing would force z=0, contradicting E_m(0)=1.
This proves the exact length-two property without Schur irreducibility.

For a single nonzero monomial f=a x^r, 0<=r<=2,
Lambda(f^n)=a^n(rn)! is nonzero for every n>=1.
These cases and the dense case exhaust all nonzero supports of degree<=2,
proving Theorem 1 with the actual support size in every case.

This completes only the univariate quadratic theorem. It does not prove
the general six-monomial binary quadratic theorem, standard SFC(2) in
two variables, or full SFC(3). The already proved reversible adjunction
bridge can carry the dense part to its stated four-term binary and
five-term ternary classes; it does not justify an unrestricted extension.
