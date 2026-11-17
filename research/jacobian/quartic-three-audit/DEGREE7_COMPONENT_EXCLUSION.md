# Degree-seven cubic-component exclusion

Proof audit dated 15 November 2026. This extends GLDC9 by new field,
fibre-gluing and low-degree power-class arguments. It is not an extrapolation
of the degree-six bookkeeping. The universal proof is written, with exact
symbolic controls; it has no formalization or identified external referee.

Claim: no P,Q in C[x,y,z] of total degree at most 7 have
Jac(P,Q,r) = kappa != 0 for r=x+x^2*y. Consequently no three-dimensional
Keller counterexample of degree at most 7 has an output combination of
degree at most 3. Relative to GLDC9 this closes 733,743,753,763. This is an
exclusion, not a new witness. The global 3D degree interval stays [5,7];
the cubic-component subclass becomes [8,9], with its old degree-nine control.

## Retained setup

Use B=C[x,y], h=1+xy, r=xh, u=1/x, K=C(r), F=K(u)=C(x,y),
and D=-x^2*d_x+(1+2xy)*d_y=d_u. Write
[P,Q]=DP*Q_z-P_z*DQ. Its sign can be absorbed in kappa.
GLDC8/9 supply the polynomial kernel B intersect K=C[r], squarefree fibres,
no polynomial slice, monic constant-leading obstruction, independent-of-z
obstruction, all-degree affine-pair obstruction, affine fibre-factor lemma,
and no-fibre affine-leading obstruction. All degree-increasing operations
below are followed by degree-free arguments; no old total-degree bound is
silently reused.

## Power classes and an invariant-leading restriction

For a polynomial a in K*F^k, UFD valuations give a=H(r)*s^k up to a scalar,
where H is a polynomial with exponents in {0,...,k-1} and s in B. Indeed,
non-fibre prime valuations are multiples of k; within each squarefree fibre,
all component valuations have the same residue modulo k. This includes the
two factors x,h at r=0. The selected H divides a, and the quotient is a
polynomial kth power. In particular:

* degree(a)<k and degree(a)<=3 imply a is invariant;
* degree(a)<=k+1 for k=3,4 gives either invariant a or a scalar times an
  affine-linear kth power (the latter only if its degree fits);
* the square-class list at degree<=5 is s^2, degree(s)<=2, or
  (r-c)*ell^2 with ell affine-linear, including constant ell.

If p^n/q^m is in K, set g=gcd(m,n). K is algebraically closed in K(u),
so p^(n/g)/q^(m/g) is in K. Bezout on the coprime reduced exponents implies
p is in K*F^(m/g) and q in K*F^(n/g).

A further restriction: a polynomial component with invariant leading
coefficient p(r) can occur in a constant-bracket pair only if r divides p.
For S=P^(-1/m), the formal residue identity from GLDC8 gives, after multiplying
by the invariant p^(1/m), a G in B[1/p] with DG in C*. Clear a finite power of
p to obtain g0 in B with Dg0=c*p^N, c!=0. A polynomial solution of
Dg0=L(r) requires r|L. To see this, write g0=L(r)/x+E(r). The elementary
intersection lemma in the next paragraph gives E polynomial and r|L.
Thus p(0)=0. In particular every nonconstant affine-in-r leading coefficient
is a scalar multiple of r.

Intersection lemma: A(r)/x+B(r) in B implies A,B in C[r] and r|A.
At a denominator root c, on r=c with c!=0, x is a free nonconstant invertible
coordinate. The leading pole coefficient A0/x+B0 cannot vanish identically
unless both constants vanish. At c=0 use the component h=0, where x is
again free and invertible, giving the same conclusion. After removing all
finite denominator roots, A,B are polynomials; reduction modulo x then
forces A(0)=0. This proof also justifies its use for L above.

## Rational monic pairs of z-degrees (3,2) and (4,3)

The GLDC9 (3,2) computation shows that monic cubic/quadratic rational pairs
in F[v] cannot have Jacobian J in K*. For J=kappa*u/H(r), the same equations
force b^2=-3*kappa*u^2/(2H)+E(r). Rational squareness forces E=0 and
-3*kappa/(2H) to be a square in K. A nonconstant affine H has a simple
zero, so this is impossible. These are field necessities, not polynomial
target transformations.

For a monic (4,3) pair, depress the cubic, and subtract a K-multiple of Q
from P. Coefficient comparison forces

```
Q=v^3+A*v+B,
P=v^4+(4*A/3+e)*v^2+(4*B/3+f)*v+2*A^2/9+2*e*A/3+g,
e,f,g in K.
```

Put X=A+3e/2, Y=B+3f/4. The v coefficient of the bracket is
-4*(XY)'/3, and the constant coefficient is the derivative of
T=4*X^3/27-e*X^2/3-2*Y^2/3. Thus XY=H in K and T=J*u+L if J in K*.
If H!=0, substitute Y=H/X: the rational function
4*X^3/27-e*X^2/3-2*H^2/(3*X^2) has degree 5 as a map of P1, whereas
composition with a nonconstant rational X(u) must have degree a positive
multiple of 5, not 1. If H=0, either X=0 and T is a nonzero scalar times
Y^2, or Y=0 and T is a cubic polynomial in X. Both composition degrees
are at least 2. Constant X,Y also cannot give J!=0. This excludes (4,3).
All uses of rational-map degree can be made over an algebraic closure of K.

## The invariant-leading (4,2) case

Suppose both leading coefficients are nonconstant affine functions of r.
The invariant-leading restriction and constant output scalings make both r.
Write Q=r*z^2+a*z+b. Subtract Q^2/r from P over F[z]. A nonzero residual
cubic coefficient belongs to K and would give the impossible monic (3,2)
case with J in K*. Remove the residual quadratic by a K-multiple of Q.
Depressing Q, the residual affine pair then forces the following normal form:

```
v=z+a/(2*r),
Q=r*(v^2+alpha(r)*u+beta(r)),
P=Q^2/r+H(r)*Q+C(r)*v+E(r),
C=-kappa/(r*alpha), alpha!=0.
```

These are rational identities for globally polynomial P,Q. Polynomiality of Q
and the intersection lemma give r^2*alpha=r*A(r), r^2*beta=B(r), both polynomial,
and

```
a^2=4*r*b-4*A(r)*h-4*B(r).
```

The z^2 coefficient of P is a^2/r+2*b+H*r. Consequently T=r^2*H is polynomial
and a^2+T is divisible by r. On the two disjoint components of r=0, h takes
values 1 and 0. The last two identities force A(0)=0, so alpha is polynomial.
They also give a^2 modulo r equal to a constant. The two component values of
a are therefore constants c1,c2 with equal squares, and T(0)=-c_i^2.

The z coefficient of P is 2*a*b/r+H*a+C. Let t=ord_r(alpha)>=0.
C has pole order t+1; the other two terms have orders at most 1 and 2.
For t>=2 no cancellation is possible. If t=0, any nonzero c_i gives an
uncancelled double pole from H*a. Thus c1=c2=0, r|a and r|T, after which
the other two terms are regular and the simple pole of C is impossible.
If t=1, cancellation of the double poles requires
T(0)*c_i=kappa/alpha1(0), where alpha=r*alpha1. Thus c1=c2=c!=0,
alpha1(0)!=0. A polynomial source z-shift makes a=c, preserving both
leading coefficients and the bracket. Then b=r*alpha1*h+gamma(r), with
polynomial gamma by the kernel/intersection lemma. Subtract gamma from Q.
The z^2 coefficient of P forces H=-c^2/r^2+L(r)/r with L polynomial.
Its constant coefficient, after subtracting polynomial terms, is

```
-c^2*alpha1(r)*h/r + I(r) = -c^2*alpha1(r)/x+I(r).
```

The intersection lemma says r|alpha1, contradicting alpha1(0)!=0.
This excludes the last invariant-leading case without bounds on a,b.

## Complete degree-seven split

Order m=deg_z P>=n=deg_z Q>0. Both leading coefficients p,q are nonconstant,
so m,n<=6, degree(p)<=7-m, degree(q)<=7-n. The highest relation is
p^n/q^m in K. All 21 integer pairs are covered as follows.

Equal degrees: for m>=5 both coefficient degrees are <=2, so their ratio
is constant and a constant target shear lowers the pair. For m=4,
a nonconstant ratio forces both coefficients affine in r; a constant target
combination gives constant leading coefficient or lowers degree. For m=3,
a nonconstant ratio has p=eta*A(r),q=eta*B(r), degree(eta)<=1 and independent
affine A,B. If eta is constant use the preceding argument. Otherwise normalize
the two leading coefficients to eta,r*eta, then subtract r times the first
component from the second. Its z-degree is <=2, while the first is cubic with
linear leading eta. A positive remaining second degree gives eta in K*F^3,
impossible by the power-class lemma; zero degree is already forbidden.
For m=2 the same argument has degree(eta)<=2. After the r shear the second
component is affine in z. The affine fibre lemma and degree(eta)<=2 force
eta/q1^2 to be constant: a numerator fibre is too large and a denominator
fibre would also divide eta. Subtract a constant multiple of the second
component squared to reach the all-degree affine-pair obstruction. Equal
m=1 is that obstruction directly.

Unequal n>=3: except (4,3), the power-class exponent m/gcd(m,n) exceeds
degree(p)<=7-m, which is <=2. Its only low-degree power-class possibility
is constant p, forbidden. For (4,3), the degree<=3 coefficient p is in
K*F^4 and hence invariant; q is invariant too. Normalize to monic over K
and use the new rational (4,3) obstruction.

Unequal n=2: q is square up to K by GLDC9's quadratic-z lemma.
For m>=4, p is in K*F^m. If m>=5 its degree<=2 makes it constant, impossible.
If m=4, its degree<=3 makes it affine invariant; q is then affine invariant,
and the new (4,2) obstruction applies. If m=3, degree(p)<=4 gives either
invariant affine p (and hence invariant affine q), excluded by monic (3,2),
or p=const*ell^3. Then q=H(r)*ell^2 with H polynomial, degree(H)<=1:
a denominator would require a complete fibre in ell^2, impossible.
If H constant use GLDC9's affine-linear residue and x^3/x^2 regularity
obstructions. Otherwise the same residue argument forces ell proportional
to x; with v=xz the monic Jacobian is const*u/H(r), impossible by the
new rational-square condition. This explicitly covers the new degree-five
square-class coefficient, rather than ignoring it.

Unequal n=1: for m>=5, p is in K*F^m and has degree<=2, impossible.
For m=3,4 the power-class lemma gives either p=const*ell^m or invariant
affine p. In the first case q=h0(r)*ell, with no denominator, and a nonconstant
h0 would force a complete fibre into ell^m by the affine fibre lemma,
impossible. Thus q has no fibre and the all-degree no-fibre lemma applies.
In the invariant case q is invariant of degree<=6. The fibre lemma forces
q=const*(r-c)^k, k=1 or2, with p=const*(r-c).
If k>=2, the next coefficient relation q*D(p_(m-1))=m*p*Db forces Db
vanish on the fibre, contradicting the constant bracket there.
If k=1 and m>=3, p_(m-1)=const*m*b+E(r); specialization forces p_(m-1)
vanish on the fibre, again Db=0 there. For m=2 the same k=1 branch is
GLDC8's polynomial-slice identity; k=2 is the preceding divisibility argument.
Finally m=2 non-invariant p has degree<=5 and is either s^2, degree(s)<=2,
or (r-c)*ell^2 with nonconstant linear ell. In the square case q=h0(r)*s;
its denominator is impossible and any numerator fibre contradicts the
fibre lemma and degree(s)<=2. A constant shear reaches affine pairs.
In the latter case q=h0(r)*ell, and polynomiality and degree<=6 make h0
at most affine. If h0 constant, use the no-fibre lemma. Otherwise the
fibre lemma forces q=const*(r-c)*ell. Normalize the two constants to get
p=f*ell^2,q=f*ell, f=r-c. The z coefficient gives
D(a/ell-2*b)=0, so a=ell*(2*b+E(r)), with E polynomial since a,ell are
polynomial and no fibre fits in ell. The constant bracket is then divisible
by ell, impossible. Constant ell is the invariant case already covered.

## Consequence and exact boundary

The degree-independent cubic scalar classification reviewed in GLDC7/8
uses only noncriticality and the exact residue two-form. It gives either
r above up to affine equivalence, or a coordinate completion with inverse
degree<=5. In the latter case the plane-fibre degree is <=35, within Moh100.
GLDC5 handles smaller scalar degrees. Thus all four remaining degree-seven
profiles with minimum entry3 are excluded, leaving12 whole profiles:
544,554,555; 644,654,655,664,665,666; 744,754,755.
No three-dimensional degree-five or-six or-seven candidate can have minimum
invariant output degree<=3. The upper witness764 and all positive sparsity
and higher-dimensional construction bounds are unchanged.

The cubic-component degree bound is [8,9]. This proof does not exclude
degree8: e.g. the low-degree power classes and the equal-degree residual
leading coefficients acquire additional possibilities. Neither absolute
3D degree nor term or special-form dimension optima are claimed.


## Replay and trust boundary

Run `.venv/bin/python research/jacobian/quartic-three-audit/verify_degree7.py`.
The [receipt](degree7-validation.json) checks exact monic and gluing identities,
all21 integer degree cases and the62/12 partition against the prior58/16
receipt. It pins the retained GLDC7/8/9 premises and the degree9 positive
control. It does not independently establish the universal UFD, rational-map,
fibre-specialization or geometric classification arguments. Those are the
written proof above and the explicitly retained source-audit premises.
Characteristic-zero coefficients and any hypothetical collision descend to
a finitely generated field embeddable in C, as in the preceding exclusions.
