# Sharp degree nine for a cubic output combination in dimension three

Dated 15 November 2026. There is no characteristic-zero three-variable
Keller counterexample of ordinary degree <=8 with a nonzero constant output
combination of degree <=3. The retained degree-nine control has a cubic
component, determinant1 and a rational collision, so the minimum in this
subclass is exactly9. This does not settle the global 3D minimum, still[5,7],
term minima, or higher-dimensional special-form dimension minima.

The new universal arguments below are written proofs. The exact checker
reconstructs the algebraic identities and the positive determinant/collision;
it does not formalize or independently referee the universal proof. The
GLDC7/8 degree-independent cubic scalar classification and KPAIR1/Moh100
remain explicit retained premises. No novelty priority is asserted.

## Setup and retained results

Let B=C[x,y], h=1+xy, r=xh, u=1/x, K=C(r), F=K(u),
D=-x^2*d_x+(1+2xy)*d_y=d_u. Use [P,Q]=DP*Q_z-P_z*DQ=kappa!=0.
GLDC8-10 give squarefree fibres, B intersect K=C[r], no polynomial slice,
constant-leading and z-independent obstructions, the all-degree affine-pair
obstruction, affine fibre-factor and no-fibre lemmas, the UFD power-class
classification a in K*F^k => a=H(r)*s^k with H kth-power-free polynomial,
and the invariant-leading restriction r|p(r).

The intersection lemma has a local version: if A(r)/x+B(r) is regular at
both prime components x=0 and h=0, then A,B are regular at r=0 and A(0)=0.
At h=0, x is a free invertible coordinate, so principal parts of A/x+B
cannot cancel unless both vanish; then inspect x=0. More generally a
combination A(r)+x^j*B(r), j>0, regular at these primes forces both A,B
regular at0, by the h=0 component. Globally, regularity on every fibre
forces A,B polynomial. All local arguments allow invariant denominators
that are units at0. Thus the GLDC10 invariant-leading(4,2) gluing proof
also applies when both leading coefficients have order1 at0, after division
by their invariant units. It uses only these local intersection statements.

Prime leading-class lemma: if a component has prime positive z-degree l,
its leading coefficient q belongs to K*F^l. If a mate has degree divisible
by l, subtract its K-leading multiple of the corresponding power of this
component over F[z]. Iterate. A nonzero residual degree prime to l gives
the assertion by the highest bracket relation; residual degree0 cannot
have constant bracket. This extends GLDC9's l=2 argument to l=3,5,7.
These are field necessities, not global polynomial target operations.

## New rational monic obstructions

All rational coefficients here belong to K(u), and prime denotes d/du.

For (5,2), depress Q=v^2+g, and subtract invariant multiples of Q^2,Q,1
from P as allowed by the coefficient equations. The remaining odd part is
v^5+(5g/2+b)*v^3+(15g^2/8+3b*g/2+c)*v, b,c in K. Hence
J=-(15g^2/8+3b*g/2+c)*g', and a cubic polynomial in g equals -J*u+constant.
For J in K*, rational-map degree3 cannot compose to degree1. No such pair.

For (5,3), depress Q=v^3+A*v+B and subtract a K-multiple of Q. The higher
coefficient equations give

```
P=v^5+c*v^4+5*A*v^3/3+(5*B/3+4*c*A/3+f)*v^2
  +(5*A^2/9+4*c*B/3+g)*v+10*A*B/9+2*c*A^2/9+2*f*A/3+j,
c,f,g,j in K.
I=5*A^3/27-4*c*A*B/3-g*A-5*B^2/3-2*f*B,
T=5*A^2*B/9+4*c*A^3/27+f*A^2/3-2*c*B^2/3-g*B.
```

The remaining coefficients are I'=0 and T'=J. After
Y=B+2*c*A/5+3*f/5, the curve I=constant becomes Y^2 equal a cubic in A
with leading coefficient1/9. This curve is irreducible over an algebraic
closure of K and has one smooth point at infinity, with pole orders2,3
for A,Y. T has unique highest-pole term5*A^2*Y/9, of order7; every other
term has order<=6. On the normalization, T therefore has degree7 as a map
to P1, including all singular specializations of the finite cubic.
If A,B are nonconstant rational functions of u, their map P1_u to this
normalization is finite onto it, and T(u) has degree a positive multiple
of7. It cannot equal J*u+constant with J in K*. Constants cannot give J!=0.

For (5,4), depress Q=v^4+A*v^2+B*v+C and subtract a K-multiple of Q.
The higher coefficients force

```
P=v^5+(5*A/4+b)*v^3+(5*B/4+c)*v^2
  +(5*C/4+5*A^2/32+3*b*A/4+d)*v
  +5*A*B/16+3*b*B/4+c*A/2+e,
b,c,d,e in K,
I1=5*A^3/32+3*b*A^2/8-5*A*C/4-d*A-5*B^2/8-2*c*B-3*b*C,
I2=15*A^2*B/32+c*A^2/2+3*b*A*B/4-5*B*C/4-d*B-2*c*C.
```

The remaining equations are I1'=I2'=0 and
omega=B*d(5*A*B/16+3*b*B/4+c*A/2)
 -(5*C/4+5*A^2/32+3*b*A/4+d)*dC=J*du.
Work over an algebraic closure of K. If rational A,B,C are not all constant,
at least one has a pole on P1_u. At any such point:

* A finite and B infinite is impossible: I2 forces C bounded, then the B^2
  term of I1 cannot cancel.
* A,B finite and C of pole order k>0 gives omega pole order2k+1>=3,
  from -5*C*dC/4; the other differential is regular.
* A,B both infinite: I2 gives C~3*A^2/8; I1 then gives B^2~-A^3/2.
  Set A~a*t^(-2k), B~b0*t^(-3k), b0^2=-a^3/2. The leading coefficient of
  omega is55*k*a^4/32 at order t^(-8k-1), nonzero.
* A has pole order k>0 and B finite: I1 gives C~A^2/8. The -dC term gives
  omega leading5*k*a^4/64 at order t^(-4k-1), while B*d(...) has order
  at most k+1. Again the pole has order>=5.

But J*du, J in K*, has no finite pole and only a double pole at infinity.
Every possible pole is contradictory. Constants give J=0. Thus (5,4) is
also impossible. This argument includes zero values of b,c,d and every
finite degeneration, without divisions by those constants.

## Two common linear-factor obstructions

The GLDC9 residue argument for leading coefficient ell^m gives a rational
primitive of1/ell whenever ell is nonconstant affine-linear. Its complete
residue strata force ell proportional to x. Apply this before the two
new common-power cases below; constant factors are absorbed by scalings.

For original leading coefficients x^4,x^3 and z-degrees(4,3), put v=xz
and depress the cubic. The monic Jacobian is J=kappa*u. The GLDC10 equations
XY=H and T'=J now give T=kappa*u^2/2+L. The degree5 branch H!=0 and the
cubic branch Y=0 cannot compose to degree2. Thus X=0 and Y=lambda*u,
lambda in C*, after the rational-square test. In particular the depressed
cubic has A in K and B=lambda*u+E(r). Before depression write
Q=x^3*z^3+q2*z^2+q1*z+q0. The invariant A condition is
x^3*q1-q2^2/3=x^4*A(r). Polynomiality makes A polynomial (its denominator
would need a full fibre in x^4). Valuation at x=0 forces x^2|q2 and x|q1.
Then the depressed constant B=q0-q2*q1/(3*x^3)+2*q2^3/(27*x^6) is polynomial.
This contradicts B=lambda/x+E(r) by the intersection lemma.

For leading x^4,H(r)*x^2 of degrees(4,2), H is nonzero constant or affine.
Put v=xz+q1/(2*H*x), so Q/H=v^2+g. The monic quartic equations give

```
P=(Q/H)^2+L(r)*(Q/H)+a(r)*(v^3+3*g*v/2)+C(r)*v+E(r),
(3*a*g/2+C)*g'=-kappa*u/H.
```

If a!=0, rational squareness after integrating forces

g=alpha(r)*u+beta(r), a*alpha^2=-2*kappa/(3*H),
and 3*a*g/2+C=lambda(r)*u, lambda=-kappa/(H*alpha).
Let U=H^2*alpha,V=H^2*beta. Polynomiality of Q gives U,V polynomial,
x|q1, q1=x*t, and r|U, via
H*q0-t^2/4=U/x+V. The z^3 coefficient of P is
x^3*(2*t/H+a), a=-2*kappa*H^3/(3*U^2).
If H(0)!=0 its pole at0 cannot cancel. Thus H is proportional tor;
normalize H=r. If k=ord_r U, k>=3 gives an uncancelled pole in a.
If k=1, the z^3 coefficient forces h|t, not necessarily r|t.
Write t=h*s and U=r*A1 with A1(0)!=0. From the Q identity, V=r*beta and
q0=h*s^2/(4*x)+A1/x+beta. The z^2 coefficient of P is
x^2*(t^2/r^2+2*q0/r+L+3*a*t/(2*r)). At h=0 all terms except L and
2*q0/r are regular, so L has at most a simple pole and its residue would
have to cancel 2*A1(0)/x+2*beta(0). An invariant residue cannot cancel the
nonconstant1/x. Thus k=1 is impossible.
Hence k=2; alpha=U/r^2 is polynomial and a unit at0. Any other zero of alpha
also gives an uncancelled pole in the z^3 coefficient, so alpha is a nonzero
complex constant. Put c=kappa/(3*alpha^2). The z^3 coefficient forces
h|(t-c). Define w=xz+(t-c)/(2*r), which is rational but regular along h=0.
Since V(0)=-c^2/4, gamma=(V+c^2/4)/r is polynomial. Subtract gamma from Q:

```
Q=r*w^2+c*w+alpha*h,
v=w+c/(2*r),
P=(Q/r)^2+L*(Q/r)-2*c*v^3/r-3*alpha*c*u*v/r+E.
```

No assertion that w is a polynomial coordinate is needed. Its constant
z coefficient is regular along h=0. There x is a unit, and the z^2
coefficient forces L=2*c^2/r^2+b(r), with b regular at0. The z coefficient,
after terms regular at h=0, has
c^3*x/(2*r^3)+c*b(r)*x/r-alpha*c/r. Its triple pole cannot cancel.
This excludes a!=0, including opposite signs on the two fibre components.

If a=0, g=A(r)*u^2+B(r), A!=0, C=-kappa/(2*H*A).
If H(0)!=0, work locally at0 and divide Q by its invariant unit H.
Polynomiality of Q implies A,B regular at0 by the local intersection
statements. The z^2 and z coefficients of P successively force L and C
regular at0. Since A*C is a unit, both A,C are units. The constant
coefficient at h=0 makes E regular at0; its value at x=0 then forces x|q1.
The identity defining g implies ord_r(A)>=2, contradicting that A is a unit.
If H=r, the z^3 coefficient initially gives only h|q1, so write q1=h*t.
The Q identity is
x^4*q0=r*A*x^2+r*B*x^4+r*t^2/4. Its regularity on every fibre, using h=0
at r0, makes U=r*A and V=r*B polynomial. If U(0)!=0, its x-adic order2
on the left cannot equal the odd order1+2*ord_x(t) on the right.
Thus r|U, so A is polynomial; the same equation then forces x|t.
Write t=x*s, giving q1=r*s and q0=r*(A+s^2/4)/x^2+V.
Subtract V from Q. The z^2 coefficient forces the adjusted L polynomial.
The z coefficient is then 2*s*(A+s^2/4)/x^2+L*s+C*x.
At h=0 the first two terms are regular, while C=-kappa/(2*r*A)
has a pole because A is polynomial. Contradiction.
Thus every new (4,2) common-power branch is excluded.

## The second invariant-leading (4,2) branch

For invariant leading p of degree<=4 and q of degree<=6, the invariant-leading
restriction gives p=const*r and q=r*(a*r+b). If b!=0, locally divide out
units and apply the GLDC10 order(1,1) proof. It remains p=r,q=r^2.
The rational normal form is
Q=r^2*z^2+q1*z+b, P=Q^2/r^3+H(r)*Q+C(r)*v+E(r),
v=z+q1/(2*r^2), g=alpha*u+beta. A residual cubic is excluded by monic(3,2).
The z^3 coefficient is2*q1/r, so q1=r*a. Polynomiality of Q gives
b-a^2/4=A(r)*h+B(r) with A,B polynomial. The z^2 coefficient gives
T=r^3*H polynomial and a^2+2*b+T divisible by r. On each component of r0,
these two relations make a^2, then a and b, constants. Thus DQ and Q_z both
vanish modulo r, contradicting the constant bracket. This closes order(1,2).

## Affine mates and equal-degree residuals

For a quadratic-z P with invariant leading p(r) and affine Q=q(r)z+b,
the fibre lemma gives every root of q as a root of p. If a root has higher
multiplicity in q than p, q*Da=2*p*Db forces Db vanish on that fibre,
contradicting the bracket there. Thus q divides p, and a=2*(p/q)*b+E(r)
with E polynomial. Then kappa=D(q*c-(p/q)*b^2-E*b), a forbidden polynomial
slice. This excludes every invariant-leading quadratic/affine pair.

At total degree8 an affine mate with m=2 has leading p of degree<=6.
The power-class list is s^2 (degree(s)<=3), (r-c)*ell^2 with linear ell,
or invariant p. If non-invariant s contains a complete fibre, degree3
would make it invariant, impossible. The no-fibre lemma and fibre-factor
lemma therefore close the square case. In the mixed case q=h0(r)*ell,
with h0 polynomial of degree<=2 and roots only c. For h0 constant use
no-fibre; for h0=const*(r-c) use GLDC10's ell-divisibility identity.
For h0=const*(r-c)^2, the next coefficient gives
(r-c)*(ell*Da-a*Dell)=2*ell^2*Db. Choose an irreducible component of r-c
not dividing ell (one exists because a complete fibre has degree3).
Then Db and Q_z vanish there, impossible. This includes ell=x,c=0 by
choosing the h=0 component; no illicit coprimality with the whole fibre.
For m=3,4,5 the leading power-class list is unchanged in kind: ell^m or
invariant affine p. The GLDC10 affine-fibre argument applies, allowing
arbitrary invariant powers in q; m>=6 has no nonconstant leading option.

For equal z-degrees, constant leading ratios lower the pair by a constant
shear, preserving degree8. For m>=5 the nonconstant ratio case has invariant
affine coefficients, and hence cancels. For m=4 or3, write p=eta*A(r),
q=eta*B(r) with independent affine A,B and degree(eta)<=1 or2 respectively.
Normalize to eta,r*eta and subtract r times one component from the other.
The residual degree is smaller. For m=4 every positive residual degree
would make eta a kth power up to K with k>=2, impossible at degree1.
For m=3 the residual gives a cube class, impossible at degree<=2.
Zero residual degree is forbidden too. For m=2, a nonconstant ratio with
max(degree(A),degree(B))=1 gives degree(eta)<=3. After the same shear its
mate is affine; eta/q1^2 in K makes eta either a linear square or invariant
affine. The former closes by the no-fibre/affine-pair lemmas, the latter
by the invariant-leading affine-mate result just proved. If that maximum
is2, eta is constant: p,q are invariant of degree<=6. Each is divisible
by r; independent ones span r,r^2. An invariant shear leaves one quadratic
with leading r and the other affine, already excluded. Proportional leading
coefficients cancel by a constant shear. Equal degree1 is the affine-pair
obstruction. All increasing shears are followed by degree-free exclusions.

## Exhaustion and scope

Order m>=n>0. Nonconstant leading coefficients imply m,n<=7; there are28
pairs. The equal and affine-mate cases are covered above. For n=2 or3,
the prime leading-class lemma makes p a full mth power up to K. Degree(p)
<=8-m leaves: invariant leading pairs(5,2),(5,3),(4,3),(4,2), or leading
ell^4 for (4,3),(4,2), or the GLDC10 (3,2) square-class cases. Invariant
(4,3) uses GLDC10's monic43 obstruction. The new monic52/53 proofs and
invariant42/common-power arguments handle the others. For n>=4, the only
new survivor of the highest-relation power bound is(5,4) with invariant
leading coefficients, now closed by monic54. Every other unequal pair has
no nonconstant power-class leading coefficient within its degree bound.
This covers all28 pairs and proves no degree<=8 completion of r.

For an arbitrary scalar combination of degree3, retain the degree-independent
GLDC7/8 classification into this exceptional cubic or a coordinate with a
completion of inverse degree<=5. The latter gives plane degree<=40,
within Moh100. GLDC5 handles lower scalar degrees. Hypothetical collisions
descend to a finitely generated characteristic-zero field embedded in C.

The separately retained control is recomputed directly from its coefficient
list by the checker: component degrees(3,7,9), determinant1, distinct rational
points(0,1,-8),(-1,-1/3,4), common image(0,1,0). This proves the matching
upper9. It has19 nonlinear occurrences and supplies no global term minimum.
No additional strict degree<=7 whole profile is closed beyond GLDC10:
62 excluded,12 remaining; global3D degree[5,7], raw13, normalized[10,48],
ordinary cubic12/36 and cubic-linear89 remain unchanged and not globally
optimal. The full simplicity goal remains open.

Run `.venv/bin/python research/jacobian-followup/quartic-three-audit/verify_degree8.py`.
The [exact receipt](degree8-validation.json) is an algebraic and positive-witness
control, not a formalization of all fibre/UFD/pole arguments or the imported
geometric scalar classification. All such trust boundaries remain explicit.
