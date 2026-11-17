# Polynomial quintics with high contact against the complete gap21 family

14 November 2026. PF2R7B5RAT1 excludes every polynomially normalized
degree(3,5) member of the closest-quintic pencil from PF2R7B5CONTACT1.
Consequently finite intersection multiplicity is at least10, and infinity
contact41 is impossible. PF2R7B5REC1 gives exactly seven degree(3,5)
polynomial-image candidates, counted with multiplicity, for infinity
contact at least38 against each first image C. There are seven distinct
candidates on a nonempty open in each family.

The same-coordinate degree(3,5) hypothesis and the additional high-contact
condition are explicit. Neither is derived for the actual unknown second
ramified image. These candidates are actual polynomial image curves;
their compatibility with Keller covering data is not established.
Written algebraic geometry and exact identities, without formalization
or independent external review.

## 1. Fixed coordinates and all normalization freedoms

Retain the complete L and scaled K charts, parameters d and nonzero z,
from PF2R7FAM1. Let G_C be the explicit closest quintic in
[PF2R7B5CONTACT1](../jc2-r7-complement/PROOF.md).
It is monic in V^3 and starts with V^3-U^5.

For a polynomial image B of degrees(3,5) with this same leading ratio,
an affine change of its source parameter makes its normalization uniquely

    P(s)=s^3+A*s+M,
    Q(s)=s^5+b*s^4+c*s^3+e*s^2+f*s+n.                 (1)

The ratio condition is forced by contact at least38: the restriction
pole basis in PF2R7B5CONTACT1 permits no mismatch at degree30. If the
original leading coefficients are u,v, then v^3=u^5, so a unique source
scale simultaneously makes both leading coefficients1. Source translation
removes the quadratic term of P. Both target translations and the shear
coefficient c remain in(1); none are discarded.

Every pair(1) is birational onto its image, since the extension degree
divides3 and5. Its affine normalization is A1: s is integral over C[P],
and C[s] is integrally closed with the same fraction field. The implicit
equation is the signed resultant already proved in PF2K1I1, after the
explicit substitutions

    a=A, b=b, c=e, d=f-c*A,
    target_P=U-M, target_Q=V-c*(U-M)-n.              (2)

We take the sign making V^3 monic. Formula(2) retains the fixed target
coordinates. Uniqueness of(1) follows because any two normalization
parameters differ affinely; centering forces the translation zero, and
the simultaneous cubic/quintic leading conditions force scale1.
Distinct solutions(A,M) below therefore give distinct B curves.

## 2. Triangular recovery and the exact two-equation interface

Substitute(1) in G_C. The coefficients at s^14,s^13,s^12,s^11,s^10
are successively linear with constant coefficient3 in b,c,e,f,n.
Thus they uniquely recover those five coefficients polynomially from
A,M and the chosen parameter of C. In particular c=5*A/3. Write

    G_C(P(s),Q(s))=sum_(j=0)^9 r_j*s^j.

The complete recovery formulas and r9,r8 are in
[normalizations.json](normalizations.json). Construct the actual implicit
quintic F_(A,M) by(2). A direct polynomial identity gives

    G_C-F_(A,M)
       =r9*U^3+r8*U*V+k6*U^2+k5*V+k3*U+k0.         (3)

The coefficient polynomials k_j are supplied explicitly in that JSON.
The checker constructs both sides from the original parametrizations,
rather than assuming ideal membership from a printed basis.
In particular r7=3*A*r9+b*r8 identically.

The pole degrees on C of U^3,U*V,U^2,V,U,1 are respectively
18,16,12,10,6,0. Since G_C restricts to degree9, equation(3) proves:

* An actual polynomial degree(3,5) curve has infinity contact at least38
  with C if and only if it arises from a solution r9=r8=0.
* At such a solution its actual implicit equation is
  G_C-k6*U^2-k5*V-k3*U-k0.
* Contact is38 when k6 is nonzero,40 when k6=0 and k5 is nonzero.
  The remaining possibility41 is excluded in Section4.

Necessity uses the thirteen-element pole basis: contact at least38 means
finite restriction degree at most12, so all higher basis coefficients
must vanish. Conversely(2) constructs a genuine polynomial normalization
at every solution; equation(3) computes its contact. This is stronger
than recognizing only rationality or total delta.

## 3. Seven candidates at every parameter, including exceptional C

Work over QQ[d] or QQ[z], without inverting any parameter. Give A,M
weights2,3 and order monomials by(weight,A-exponent). Put

    f0=(27/5)*r9, g0=(9/5)*r8,
    h0=-(M*f0-A*g0)/9.

For each family these three polynomials have monic leading monomials

    A^3, A^2*M, M^3.

Their three S-polynomials reduce to zero over the polynomial coefficient
ring. This is checked by exact monic division. Since h0 is a polynomial
combination of f0,g0, this is a Groebner basis of exactly their ideal,
with no hidden localization or extra equation. Its leading coefficients
are constant1, so every specialization has the same basis and standard
monomials

    1,A,A^2,M,A*M,M^2,A*M^2.

The quotient of the normalization coefficient space is therefore free of
rank7 over the entire parameter ring. Multiplicities refer to this
explicit coefficient scheme, not to an asserted scheme structure on
abstract curves.
Every complex fibre has length7 and is nonempty, counting algebraic
multiplicity. This includes the nonimmersed and all exceptional first
images; z0 is not in the Laurent chart. No reducedness is asserted at
every parameter.

At d1 and separately z1, exact lexicographic elimination gives a degree7
squarefree polynomial in A and a monic linear formula for M. Both are
retained in normalizations.json and replayed. Thus each family has a
nonempty open with seven distinct choices. The finite free algebra,
not a numerical parameter sweep, supplies the uniform count.
Sections1--3, with Section4's contact restriction, prove PF2R7B5REC1.

## 4. No degree(3,5) polynomial image has contact41

Contact41 would put the actual implicit equation in

    G_C+alpha*U+beta.

Substitution in(1) first forces the same five recovery equations. Since
alpha*P+beta has degree at most3, every residual r9,r8,r7,r6,r5,r4
must vanish. Primitive integer multiples of these six residuals generate
the necessary ideal I. Frozen certificates prove:

* In L, a literal integer combination of the six generators is a nonzero
  constant, for all d.
* In K, two literal integer combinations give H1(z),H2(z). A retained
  rational Bezout identity gives a*H1+b*H2=z^7.

Every cofactor identity is replayed by direct convolution of integer
coefficients using Python Fraction arithmetic. The univariate Bezout
identity is replayed the same way. Therefore no solution exists for any
L parameter or any nonzero K parameter. Generic rational-function
denominators are not discarded: the second identity and Bezout replay
cover all their exceptional zeros. The discovery solver's unit output
alone is not evidence for this conclusion.

Together with PF2R7B5CONTACT1, the complete necessary infinity-contact
spectrum is now contained in

    {20,22,24,26,28,30,32,34,38,40},

and total finite intersection multiplicity is at least10. This is
PF2R7B5RAT1. It is not a bound on the number of distinct target values.

## 5. What is and is not resolved

The maximal41 pencil has no member with the required polynomial
normalization. The next high-contact candidates are explicitly constructed
as a finite free degree7 correspondence. Their contact40 sublocus is
k6=0; its full parameter stratification has not been certified here.
Nor are all seven choices asserted immersed, pairwise distinct at every
parameter, or compatible with the omitted-prime ledger.

An actual additional image is not proved to have degrees(3,5), or to
have contact at least38. Other contact levels, other target charts,
the optional unramified image, joint complement relations, boundary
attachment and source A2 realization remain unresolved. The individual
first-image topology theorem supplies none of these missing conditions.
No whole terminal-r7 or JC2 closure follows.
