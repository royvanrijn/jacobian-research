# Degree-six cubic-component exclusion

This proof excludes every degree<=6 completion of r=x+x^2y with a nonzero
constant Jacobian. Together with the degree-independent scalar classification
reviewed in GLDC7/8 and Moh100, it excludes every3D degree<=6 counterexample
with an output combination of degree<=3. The universal UFD/fibre/field
arguments are written; exact symbolic identities and all15 degree cases
are separately checked. No formal verification, external refereeing or
novelty priority is claimed. The original derivation draft remains retained.

Relative to the prior54/20 dense partition, the additional full profiles
are633,643,653,663, giving58 excluded and16 remaining. The3D ordinary-degree
interval remains[5,7]. A checked degree9 source-shear control gives only
an upper9 in the cubic-component subclass, now with lower7. No smaller
actual witness or global simplicity optimum follows. Dated15 November2026.

## Fields and basic lemmas

Put B=C[x,y], r=x+x^2y, u=1/x, K=C(r), F=K(u). The tangent derivation
D=-x^2*d_x+(1+2xy)*d_y is d_u. Its polynomial kernel is C[r]. Every fibre
r-c is squarefree; for c!=0 it is irreducible of degree3 (linear in y with coprime
coefficients x^2 andx-c), and r=x*(1+xy)
at c=0. The no-polynomial-slice, monic-in-z formal-residue and all-degree
affine-pair obstructions are exactly those proved in GLDC8.

Fibre-factor lemma: if Q=qz+b and P has z-degree>=2, then any complete
fibre r-c dividingq divides the leading z coefficient ofP. In each invariant
irreducible quotient of the squarefree fibre, the bracket becomes-P_z*Db.
A nonzero constant unit product in a domain polynomial ring has z-degree0.

No-fibre affine leading coefficient lemma, valid in any total degree:
if Q=qz+b andq contains no complete fibre r-c, no polynomial P can have
[P,Q]_D=constant nonzero. A leading P coefficientp satisfiesp/q^m inK.
A reduced denominator H(r) must divide q^m, hence have a complete fibre
factor, impossible. Thus p/q^m is inC[r]. Subtract that invariant polynomial
multiple of Q^m fromP, preserving the constant bracket and polynomiality.
Iterate until P is independent ofz, contradicting the no-slice obstruction.
No total-degree bound is reused after these nonlinear target shears.

## A quadratic z component must have square leading coefficient overK

If Q has z-degree2 and lcq, and P has even z-degree2k, the highest equation
gives p/q^k inK. OverF[z] subtract that K-multiple ofQ^k to lower the degree.
If an odd degree remains, p^2/q^m inK implies q is a square up toK*. If
only degree0 remains, the constant bracket is impossible. Therefore every
completion requires q inK*F^2. Rational target coefficients are allowed for
this field necessity; no claim of global polynomial target transformation
is made in this paragraph.

Polynomial square-class classification: for a nonzero polynomial a inB with deg a<=4,
a inK*F^2 iff a=H(r)*s(x,y)^2 up to a nonzero complex constant, with H
squarefree polynomial inr. This follows from UFD parity. Nonvertical factors
have even exponents. In a fibre r-c, the exponents on all irreducible
components have the same parity, because r-c is squarefree; the r=0 pair
x,1+xy is included. Choose H as the product of the odd fibres, which dividesa.
Thus at degree<=4 either a is a polynomial square with degs<=2, or a is
an affine nonconstant polynomial inr ands is constant. This parity argument covers every fibre, including the reducible zero
fibre, without choosing a generic irreducibility stratum.

## Two monic cubic/quadratic rational normal-form obstructions

OverK(u), for monic degrees3 and2 in v, depress the quadratic by a rational
v-shift. Write Q=v^2+g(u), P=v^3+a(u)v^2+b(u)v+c(u). If their Jacobian
in(u,v) is J(u), coefficient comparison gives

```
a'=0, 2b'-3g'=0, 2c'-2a*g'=0, -b*g'=J.
```

Subtract a*Q fromP. Then a is inK, c becomesconstant inK, b=(3/2)g+b0,
and (b^2)'=-3J.

1. If J is a nonzero element ofK, b^2 is a nonconstant linear polynomial
in u, which cannot be a rational square. No such pair exists.

2. If J=kappa*u with kappa inC*, then b^2=-(3/2)kappa*u^2+H(r).
A rational square of a polynomial has no finite denominator. The quadratic
has a simple root over the algebraic closure ofK unlessH=0; thus rational
squareness requiresH=0. Comparing coefficients then gives b=lambda*u
with nonzero lambda inC.
Thus g=alpha/x+G(r), where x=1/u, alpha!=0.

For the global polynomial pair with original z-leading coefficientsx^3,x^2,
set v=xz. The bracket equation becomes Jac_(u,v)=kappa*u. Original
Q=x^2 z^2+q1*z+q0 gives g=q0-q1^2/(4x^2). Hence

```
x^2*q0-q1^2/4 = alpha*x+x^2*G(r), alpha!=0.
```

Polynomiality forces G(r) polynomial: its reduced denominator composed withr
would divide x^2, impossible by degree3. Modulo x, q1 is divisible byx.
The left side is then divisible byx^2, contradicting the nonzero alpha*x.
This would rule out lcP=const*x^3,lcQ=const*x^2 in all xy degrees.

## Linear common leading factor in the (3,2) case

For p=const*ell^3,q=const*ell^2 with nonconstant affine-linear ell, the
formal residue S=P^(-1/3) is defined inF((z^-1)). Its coefficient gives a
rational g inF with Dg=constant/ell. Classify when1/ell has a rational
primitive in u. Write ell=A/u+B*(r*u^2-u)+C.
If B!=0,A!=0, the cubic denominator of1/ell=u/(Br*u^3-B*u^2+C*u+A)
has nonzero generic discriminant, with leading r^2 coefficient-27A^2B^2,
and simple nonzero roots have nonzero residues. If B!=0,A=0, the reduced
quadratic denominator Br*u^2-B*u+C has nonzero discriminantB^2-4BCr;
its simple poles likewise have nonzero residues, including C=0.
If B=0,A*C!=0, the residue of u/(A+C*u) is-A/C^2!=0.
A derivative of a rational function has zero residue at every finite pole.
Thus among nonconstant affine-linear ell, only ell proportional tox admits
a rational primitive. That remaining
case is excluded by the monic J=kappa*u obstruction above.

## Complete degree-six split

Normalize a hypothetical degree<=6 pair so m=deg_zP>=n=deg_zQ>0.
Constant leading coefficients and z-independent components are already
excluded. Nonconstant leading coefficients give m,n<=5. Let p,q be leading
coefficients, degp<=6-m,degq<=6-n. Their highest relation is p^n/q^m inK.
A root of the reduced numerator forces r-c|p, a root of the denominator
forces r-c|q, using squarefreeness of the fibres.

Equal degrees: m>=4 gives degp,degq<=2, so the ratio is constant and a
constant target shear lowers m+n. For m=3, a nonconstant ratio forces both
leading coefficients affine inr; a constant target combination gives a
constant leading coefficient or lowers degree, as in GLDC8. For m=2,
p=eta*A(r),q=eta*B(r) with deg eta<=1 and independent affine A,B if the
ratio is nonconstant. If eta constant, use the preceding affine-in-r
argument. Otherwise constant target changes give leading coefficientseta
andr*eta. Subtract r times the first component fromthe second to make the
second affine in z. The remaining quadratic-z component has linear leading
eta. Its ratio eta/q1^2 inK cannot have numerator or denominator fibre roots
(the fibre-factor lemma would force degree3|eta); hence eta is a constant
multiple of a polynomial square, impossible in degree1. Equal degree1 is
the old all-degree affine-pair obstruction. All degree-increasing shears
in this step are followed by degree-free obstructions.

Unequal n>=3: only(4,3),(5,3),(5,4). Here degp<=2 anddegq<=3. A numerator
fibre is impossible; a denominator fibre makesq affine inr and thenp
invariant, hence constant, impossible. A constant ratio gives the UFD
power p=alpha*l^(m/gcd(m,n)), whose exponent exceeds6-m. All are excluded.

Unequal n=2: by the quadratic square-class lemma, q is either s^2 with
nonconstant degs<=2, or affine inr. For q invariant the highest relation
makesp invariant; m>=4 thenforcesp constant. For m=3 both leading
coefficients are inK, and rescaling them to monic gives the impossible
J inK* cubic/quadratic case. For q=s^2, p/s^m is inK; its denominator
would divide s^m and thus s contains a complete degree3 fibre, impossible.
So p=h(r)*s^m. Degree bounds exclude m=4,5; m=3 forces degs=1 andh
constant. Apply the linear-factor residue and x^3/x^2 obstructions above.

Unequal n=1: for m>=4, degp<=2. Numerator fibres are impossible and
any denominator fibre, by the fibre lemma, would also dividep. The ratio
is constant, forcing p=const*q^m, impossible for nonconstantq at this degree.
For m=3, a nonconstant ratio forcesp,q affine inr; the fibre lemma forces
p=lambda*q. The next coefficient is p2=3lambda*b+h(r), Q=qz+b.
At a root ofq, unit-product specialization requiresp2 vanish on every
fibre component. Thenb is constant modulo each component, Db vanishes
there, and the constant bracket is impossible. A constant ratio forces
q linear and no complete fibre factor, so use the degree-free affine
leading lemma.
For m=2, p/q^2 inK makesp square-class as above. Ifp affine inr, q is
also affine inr; the fibre lemma givesp=lambda*q and the GLDC8 polynomial
slice identity contradicts the bracket. If p=s^2, q/s is inK with no
possible denominator because degs<=2. Thus q=h(r)*s. A nonconstant h
would give a complete fibre dividingp, hence dividing s, impossible.
Therefore the ratio is constant and subtracting a constant multiple of
Q^2 makes both outputs affine in z, again impossible.

This exhausts all15 positive pairs m>=n,m<=5. The case list is bookkeeping;
the preceding universal arguments close each branch. For other cubic scalar
types, the reviewed classification uses only degreeR<=3, noncriticality and
exactness of the fibre residue form. These properties do not bound the other
coordinate degrees. A coordinate completion of inverse degree<=5 gives plane
degree<=30, within Moh100, excluding collisions. Degree<=2 combinations are
already covered by GLDC5. Finite-generation embedding transfers hypothetical
characteristic-zero collisions toC. This proves the four stated full-profile
closures.

## Further boundary, not a theorem

At degree7 a quadratic-z leading coefficient may be (r-c)*ell^2 of degree5,
a square up toK but not a polynomial square. This case is absent atdegree6.
Known degree9 controls have both leading coefficients containing an
invariant factor. Do not infer degree7/8 exclusion from the degree6
argument. The degree7 case remains outside this theorem. Its new square-class
patterns require additional analysis; no theorem is inferred from a bounded
search or from extending the degree6 bookkeeping.

## Replay and proof boundary

```
.venv/bin/python research/jacobian-followup/quartic-three-audit/verify_degree6.py
```

[Exact controls and application receipt](degree6-validation.json) checks the
monic coefficient equations, rational chain rule, formal-residue identity,
all affine-linear leading-factor residue strata, square-class low-degree
parity cases and the complete integer degree split. It does not independently
prove the universal factorization, residue/logarithm or fibre arguments.
The [degree5 proof](EXCEPTIONAL_DEGREE5.md) supplies the retained no-slice and
squarefree-fibre lemmas; the [quartic audit](INVERTIBILITY_AUDIT.md) supplies
the reviewed degree-free cubic scalar classification. The known cubic-component
[positive control](../equal-negative-completion/source-shear-validation.json)
has degree9,19 terms and different weights, and is not a global improvement.
