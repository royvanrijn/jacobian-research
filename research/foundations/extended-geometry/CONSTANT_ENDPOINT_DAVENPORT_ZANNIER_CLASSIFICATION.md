# Constant endpoint profiles and marked Davenport–Zannier pairs

Internal research proof, 20 November 2026. This file makes no claim about
arbitrary height profiles or arbitrary Image phases. It characterizes the
constant endpoint subclass of the exact two-height profile.

Work over a characteristic-zero field K. Normalize

    p=t^a(C(x)+t^2 D(x))/x,
    a>=0 an integer, C,D in K[x], C(0)=-D(0)=1.

The sphere relation is t^2+2xy=1. Formally put

    X=sqrt(-D/C), X(0)=1,
    u=x X^a/C=x+O(x^2),
    W=X^(a-1)/C, Phi(u)=W(v(u))*v'(u),

where v is the formal inverse of u. Constant endpoint means Phi=1, or
equivalently W=u'. By the retained endpoint theorem it implies all-order
purity at every integer a. Conversely every pure unit-center phase with
even a has this property. For odd a, purity alone does NOT force Phi=1.
Rationality in this subclass follows directly, without importing the
general pure-phase rationality theorem. The identity Phi=1 gives
uprime/u=1/(xX). If a is even, u is already rational and this identity
makes X rational. If a is odd and X were genuinely quadratic over
K(x), its involution X->-X would send u to-u. It would fix uprime/u
and negate 1/(xX), a contradiction in characteristic zero. Hence u,X
are rational for every constant-endpoint profile, including a0.

## 1. Elementary classification over K

**Theorem CE-DZ.** Apart from C=1,D=-1, the following constructions are
mutually inverse over K.

(I) A normalized polynomial pair C,D of constant endpoint type at height a.

(II) Coprime monic polynomials A,B of the same degree N>=1, with a marked
simple zero A(0)=0, B(0)!=0, satisfying the following conditions. Write p0
and q0 for their respective numbers of distinct geometric roots, and put
r=p0+q0-1. Then

    deg(A-B)=N-r,                                    (DZ1)

and every A-root other than the marked zero has multiplicity at most a.
The marked simple zero is part of the data. In particular r>=1 and N>=r.

The map from(I) to(II) is obtained by writing u/u(infinity)=A/B in reduced
monic form. The reverse map is

    kappa=(A/B)'(0)!=0,
    u=(A/B)/kappa,
    P=normalized product of all roots of AB other than 0,
    H=x P^a/u,
    C=H, D=-H P^2.                                  (DZ2)

Here “normalized product” means the unique squarefree polynomial with
these roots and constant term1. It belongs to K[x], even when the roots
do not individually belong to K. This is a bijection of actual
normalized pairs, not merely of equivalence classes.

**Proof.** Write X=P/J in reduced form, with P(0)=J(0)=1. Polynomiality
forces C=H J^2,D=-H P^2, H(0)=1. Thus

    u=x P^a/(H J^(a+2)),
    W=P^(a-1)/(H J^(a+1)).                           (DZ3)

At a root of J of multiplicity j, with H-multiplicity h, u has pole order
h+(a+2)j. Therefore u' has pole order h+(a+2)j+1, while W has pole order
h+(a+1)j. These are different. Hence W=u' forces J=1, including a=0.
If P=1, the identity (x/H)'=1/H gives H'=0 and the excluded trivial pair.
Hence r=deg P>=1. Let h=deg H and delta=1+a*r-h. If delta!=0, the growth
orders of u' and W at infinity are delta-1 and delta-r-1, respectively,
which is impossible. Thus h=a*r+1 and u has a finite nonzero value c at
infinity. Dividing W=u' by u yields

    u'/u=1/(xP).                                    (DZ4)

Every root alpha of P is simple. Indeed, logarithmic derivatives of
rational functions have only simple poles, whereas 1/(xP) has a pole
with the root multiplicity of P. Define

    n_alpha=-1/(alpha P'(alpha)).                    (DZ5)

Comparing residues in

    H'/H=(P-1)/(xP)+a P'/P

shows that the H-multiplicity at alpha equals a+n_alpha. Thus
n_alpha is a nonzero integer, n_alpha>=-a, and H has no other roots.
The residue at infinity in(DZ4), or h=a*r+1, gives

    sum n_alpha=1,
    u=x product_alpha (1-x/alpha)^(-n_alpha).         (DZ6)

The zeros of u are 0, of multiplicity1, and roots with n_alpha<0,
of multiplicity -n_alpha<=a. Its poles are roots with n_alpha>0,
of multiplicity n_alpha. Consequently its reduced numerator and
denominator have the same degree N and a total of r+1 distinct roots.
At infinity, 1/(xP) has a zero of order r+1. Equation(DZ4), with c!=0,
therefore shows that u-c vanishes to EXACT order r. Writing u/c=A/B
with A,B monic gives(DZ1), with N>=r.

Conversely suppose(II). Put f=A/B and u=f/f'(0). The logarithmic
derivative u'/u has a simple pole, with nonzero integer residue, at
every root of AB, and no other finite poles. Its reduced denominator
is thus xP. Equation(DZ1) says that f-1 has exact vanishing order r
at infinity; therefore u'/u vanishes there to exact order r+1.
Since deg(xP)=r+1, its numerator is a nonzero constant. The residue1
at 0 makes that constant1, proving(DZ4).

The rational function H=xP^a/u is regular at 0 and has value1 there.
At each unmarked A-root of multiplicity m its order is a-m>=0; at
each B-root of multiplicity n its order is a+n>=0. It has no other
finite poles, so H belongs to K[x]. With C,D as in(DZ2), X=P is the
selected square root, u=xP^a/H, and(DZ4) gives u'=P^(a-1)/H=W.
This includes a=0 without treating P^(-1) as a polynomial. The two
constructions recover u,P,H and therefore C,D and A,B exactly. QED.

## 2. Exact degree and the exceptional height zero

Every nontrivial constant endpoint profile has

    deg P=r>=1, deg H=a*r+1,
    deg C=a*r+1, deg D=(a+2)*r+1,
    least sphere degree=(a+2)*(r+1).                 (DZ7)

The last equality uses the retained exact all-height sphere-degree
formula max(deg C+a-1,deg D+a+1). Both factors in(DZ7) are at least2.
Thus no nontrivial constant endpoint profile has prime least degree,
including odd a. This conclusion applies to all pure unit-center
profiles when a is even; it does not apply to all odd-a pure profiles.

For a=0 the numerator A can have no unmarked root. Its marked root
is simple, so N=1 and r=1. Equivalently all n_alpha are positive
integers with sum1. Thus the only nontrivial normalized height-zero
profiles are C=1+lambda*x,D=-(1+lambda*x)^3, lambda in K*.
Their least degree is4. This is a purity classification, not only
a minimum-degree result.

For every integer a>=1, EVERY integer r>=1 is realized over K.
Take lambda in K*, U=1+lambda*x and

    S_r=(U^r-1)/(U-1),
    P=U*S_r/r,
    H=r^(1-a)*U^(a+r)*S_r^(a-1).                    (DZ8)

These are polynomials with constant term1. Direct substitution gives

    u=(1-U^(-r))/(r*lambda),
    u'=U^(-r-1)=P^(a-1)/H.

This is the star DZ pair with passport ((1,...,1),(r)). It is a
classical rational map, and is NOT claimed as a new witness family.
For a=0,r>1, formula(DZ8) has a negative power of S_r and does not
give polynomial H, in agreement with the complete height-zero case.
In particular for every EVEN a>=2 the nontrivial pure unit-center
least-degree spectrum is exactly

    {(a+2)*(r+1): r=1,2,3,...}.                      (DZ9)

For odd a this is the exact degree spectrum only of the Phi=1 subclass.

**Corollary CE-UnsafeSpectrum.** Every nontrivial constant endpoint
profile is unsafe, and the SAME marker x, hence the original marker Z,
is nonzero at infinitely many EVEN mixed orders. Consequently(DZ9)
is the exact UNSAFE degree spectrum in this subclass for every a>=1.
For every even a>=2 it is the complete unsafe-degree spectrum of the
entire two-height class, including zero centers. For a=0 the complete
unsafe-degree spectrum is {4}. Zero-center pure phases are all safe,
by the retained all-height zero-center theorem.

Here is a direct proof of the assertion about x; the complete retained
two-marker safety classification is not needed. Work over k=Kbar and
use a second branch y=v(-u) in the finite field k(x,y), with x=v(u).
Suppose the even series v(u)+v(-u) were a polynomial Q(u^2). Write
u=A(x)/(kappa B(x)) as above. In the discrete valuation ring

    O=k[1/x]_(1/x),

the function u is regular with value c=1/kappa!=0 at infinity. The
equation

    A(y)+kappa*u(x)*B(y)=0

has coefficients in O and leading coefficient 1+kappa*u(x), whose
residue is2. Dividing by this unit makes a monic equation, so y is
integral over O. But Q(u(x)^2) belongs to O. The identity
x=Q(u(x)^2)-y would then make x integral over O. Since x belongs to
the fraction field of the integrally closed ring O, this would put
x in O, a contradiction. Thus the even series is not polynomial.

The exact marked endpoint identity is

    Lambda(x*p^m)=kappa_(a,m)[u^m]v(u),
    kappa_(a,m)=integral_0^1 t^(am)(1-t^2)^m dt !=0

whenever am is even. Every even m satisfies this condition. Since the
even coefficients are not eventually zero, the claim follows. For a
standard balanced lift F of least degree d, the retained normalization
E2(Z F^m)=(dm+2)! Lambda(x*p^m) transfers the same nonvanishing to Z.
There is no claim that every even order is nonzero for every profile.

## 3. Closest classical statements and their exact use

The pair in(II) is exactly a Davenport–Zannier pair with a marked
simple numerator root and the stated multiplicity cap. In the
notation of Pakovich–Zvonkin, its two multiplicity partitions have
lengths p0,q0 and sum N, and(DZ1) is their extremal equality
deg(A-B)=N+1-(p0+q0). The common gcd of all parts is1 because one
part is1. Hence their standard hypotheses are satisfied.

Primary sources actually inspected:

* Pakovich–Zvonkin, *Minimum Degree of the Difference of Two Polynomials
  over Q, and Weighted Plane Trees*, arXiv:1306.4141, Theorem2.24
  (PDF p10), Theorem3.1 (PDF p11), and Proposition2.18 (PDF p9).
  [Actual PDF](https://arxiv.org/pdf/1306.4141v1).
  Theorem2.24 identifies extremal pairs with weighted plane trees;
  Theorem3.1 says a passport is realizable iff
  p0+q0<=N/gcd(parts)+1. Proposition2.18 supplies a sufficient rational
  field-of-definition condition for a unique tree; it is not a
  necessary rationality condition.
* Pakovich–Zvonkin, *Minimum Degree ... Part II: Davenport–Zannier
  pairs*, arXiv:1509.07973v2, Assumption1.1, Definition1.2 and
  Proposition1.4. [Actual PDF](https://arxiv.org/pdf/1509.07973v2).
  Proposition1.4 gives a bijection of DZ equivalence classes and
  weighted bicolored plane trees with the same passport.
* Lando–Zvonkin, *Graphs on Surfaces and Their Applications* (2004),
  Section2.6.2, equations(2.16)--(2.17), Example2.6.5, printed
  pp149--150. [Actual book](https://www-users.cse.umn.edu/~reiner/Classes/Math8680Fall2014Papers/LandoZvonkin.pdf).
  This already constructs rational Belyi maps from products whose
  logarithmic derivative has constant numerator, using integral
  residues. The mechanism in(DZ4)--(DZ6) has this classical antecedent.

Thus over C the constant endpoint profiles correspond to weighted bicolored plane
trees with a marked black leaf of weight1, other black weighted
degrees at most a, and arbitrary white weighted degrees. Equivalence
of marked pairs means simultaneous affine change fixing the marked
root, and common scaling; normalizing u'(0)=1 removes target scaling
but leaves the allowed source scaling in the profile family.

Theorem3.1 gives existence over Qbar for every such
passport satisfying N>=r. It does NOT say that all pairs, or even
some pair of every such passport, are defined over Q. Over a fixed K,
Theorem CE-DZ uses a K-defined pair WITH a K-rational marked simple
root. A K-defined unmarked equivalence class alone does not supply
that marked point. No field-of-definition obstacle is hidden in(DZ2).

This caveat is real even for a=2,r=3. Put eta^2=-2 and

    U=1+x,
    V=1+(1/2-eta/4)*x,
    L=1+(2/3+eta/6)*x.

Then

    P=U*V*L, C=U^6*V, D=-U^8*V^3*L^2,
    u=x*V*L^2/U^4, u'=L/U^5=W.                     (DZ14)

The degree is16 and the DZ passport is ((1,1,2),(4)). Complex
conjugation gives two distinct marked geometric classes. Neither
class has a real model under x-scaling: C has just two roots, with
different multiplicities6 and1, and their ratio is
1/(1/2-eta/4)=(4+2eta)/3, which is nonreal and scaling-invariant.
Yet the unmarked star is a unitree and has the Q-defined model

    A0=x^4-4x+3=(x-1)^2*(x^2+2x+3), B0=x^4.

Its two simple numerator roots are -1+-sqrt(-2), so neither is
Q-rational. This is an explicit failure of marked descent, not a
failure of the unmarked unitree theorem. In particular geometric
counts below are not automatically counts of rational families.

There is a FINITE combinatorial classification at each fixed (a,r).
If t=p0-1 counts the unmarked numerator roots, then

    r<=N=1+sum_(unmarked A-roots) multiplicity
         <=1+a*t<=1+a*(r-1).                         (DZ10)

There are only finitely many possible N and passports in this range,
and only finitely many weighted plane trees of each passport. One
marks a black leaf of weighted degree1; distinct choices are
identified only when a tree automorphism preserves the choice.
This gives ALL Phi=1 profiles over an algebraically closed field,
and hence the COMPLETE pure unit-center classification for even a,
up to the following one-parameter source scaling. After the marked
root is moved to0, the remaining affine equivalence is x->s*x,
s!=0; its action on normalized data is

    u_s(x)=u(s*x)/s,
    P_s(x)=P(s*x), H_s(x)=H(s*x),
    C_s(x)=C(s*x), D_s(x)=D(s*x).                     (DZ11)

The trivial pair C=1,D=-1 is added separately. Over K the same
description is restricted to pairs WITH the stated K-descent data;
the finite-tree list is not asserted to provide K-models automatically.
For a=1, (DZ10) forces N=r,p0=r,q0=1, so only the star passport
((1,...,1),(r)) occurs. This recovers the constant-endpoint subfamily of the retained one-pole classification
without asserting that every odd-start pure phase has Phi=1.

The passage from C to an arbitrary algebraically closed characteristic-
zero field is justified as follows, rather than by an embedding of
that entire field into C. Fix the passport, label all distinct roots,
put the marked A-root at0 and a chosen B-root at1, and require A,B to
be monic. The required coefficient vanishing in A-B, the nonzero
leading difference coefficient, and distinct-root conditions form a
finite system over Q; each nonvanishing condition is expressed with
an inverse variable. Over C it has finitely many solutions, because
the weighted-tree classes are finite and the chosen labeled roots
fix the remaining affine scaling. Let S be its finitely generated coordinate algebra over Q, including
the inverse variables. The complex algebra S tensor_Q C has finitely
many maximal ideals, by the weak Nullstellensatz and the finite
complex solution set. Stacks Project, Lemma10.61.3
[Tag0ALW](https://stacks.math.columbia.edu/tag/0ALW), then makes it
finite-dimensional over C. A Q-basis remains a basis after extension
to C, so S itself is finite-dimensional over Q. Its solution
coordinates are therefore algebraic over Q in every characteristic-
zero algebraically closed field. Base change from Qbar gives the
same finite geometric solutions; possible nilpotents create no
additional field-valued points.
Undoing the source scaling proves the claimed classification over
every algebraically closed characteristic-zero field. This argument
does not establish descent to a prescribed nonclosed field K.

## Canonical source scale and exact fixed-field descent

For a nontrivial constant-endpoint profile put

    q=lc(H)/lc(P)^a=1/u(infinity) in Kbar*.

The exact degrees are degH=ar+1 and degP=r. Under the allowed source
scale x->s*x, equation(DZ11) therefore gives q_s=s*q. Thus s=q^(-1)
produces a UNIQUE representative of each marked geometric class with
q=1. In particular the source-scaling action is free. This includes
heighta0: then q=lcH. Equivalently q is the derivative at the marked
point of the associated monic ratio A/B.

This canonical representative is functorial under every coefficient-field
automorphism: normalization uses only the rational coefficient expression
q and source substitution. Fix a characteristic-zero field K and its
algebraic closure Kbar. A marked geometric class admits a K-model if and
only if it is fixed by Gal(Kbar/K). Indeed a K-model makes the class fixed.
Conversely, if the class is fixed, every Galois conjugate of its q=1
representative is another q=1 representative of the same class, hence is
identical by uniqueness. All coefficients then lie in the fixed field K;
char0 ensures K is perfect and Kbar^Gal(Kbar/K)=K. Every K-model is
obtained from this canonical representative by a UNIQUE source scale
s in K*, namely its q parameter. No separate descent or Hilbert90
hypothesis is needed for the MARKED class.

This does not replace marked Galois invariance by invariance of the
underlying unmarked tree. A rational unmarked class may have conjugate
nonrational markings, as the explicit eta²=-2 example demonstrates.


## 4. Exact count of geometric constant endpoint classes

**Theorem CE-TreeCount.** Fix integers a>=1 and r>=1. Over an
algebraically closed characteristic-zero field, the number of
nontrivial normalized constant endpoint profiles with deg P=r,
modulo x->s*x, is

    (1/r)*binomial(a*r,r-1)
      =binomial(a*r,r)/((a-1)*r+1).                 (DZ12)

These are geometric classes: complex conjugate plane trees remain
distinct unless they are isomorphic as oriented marked plane trees.
For even a these are precisely the UNSAFE unit-center pure classes.
Zero-center pure phases are not included in this finite count.

To prove(DZ12), root a marked weighted tree at the white neighbor of
its marked black leaf, and remove the marked edge and leaf. At every
white vertex, the parent edge specifies a corner, so the remaining
incident black branches are an ordered sequence. A black branch
with k ordered white children has k+1 edges. Their positive integral
weights have sum at most a; there are binomial(a,k+1) such choices.
Let T(z) count the residual white-rooted trees by their number of
unweighted edges. The branch series and ordered-sequence construction
give, as identities of formal power series,

    B(z)=sum_(k=0)^(a-1) binomial(a,k+1)*z^(k+1)*T(z)^k,
    T(z)=1/(1-B(z)),
    T(z)=(1+z*T(z))^a.                              (DZ13)

The removed marked edge means that r corresponds to [z^(r-1)]T.
With F=z*T, equation(DZ13) becomes F=z(1+F)^a. The needed formal
Lagrange identity follows directly from residues:

    [z^r]F=Res F*z^(-r-1) dz
          =(1/r)Res F'*z^(-r) dz
          =(1/r)Res (1+s)^(a*r)*s^(-r) ds
          =(1/r)[s^(r-1)](1+s)^(a*r).

The second equality is Res d(F*z^(-r))=0; the third substitutes
s=F(z), an invertible formal coordinate. This proves(DZ12).
The marked leaf kills all plane-tree automorphisms: fixing its edge
fixes the cyclic order at its white endpoint and inductively every
edge. Thus there is no unaccounted automorphism factor. Equivalently
the scaling action on profiles is free: a stabilizer s satisfies
s^r=1 from P and s^(a*r+1)=1 from H, hence s=1.

For a=2 this is the Catalan number C_r, with initial counts1,2,5,14.
For a=1 it is1 for every r, agreeing with the star rigidity above.
For a=0, B=0,T=1, so only r=1 has one class, as already proved.

The map u/c is branched only over0,1,infinity; over1 its branch
partition is(r,1^(N-r)), with the r-fold point at the source infinity.
Only this ONE branch permutation is forced to be a single cycle.
The narrower literature terminology “single-cycle Belyi map,” when
it requires a unique ramification point over EACH of the three
branch values, cannot be applied to this entire class.

## 5. Verification scope

The attached stdlib checker verifies exact rational identities for
the star family, explicit non-star marked pairs, reconstruction in
both directions, the height-zero restriction and degree controls.
These finite checks do not replace the all-order arguments above.
The retained endpoint theorem supplies purity and marked-moment transfer;
the rationality step and elementary fixed-field iff are proved directly
above by rational-function algebra. Classical tree existence is used solely for the
tree interpretation and algebraically closed passport existence,
not for the fixed-field classification or explicit rational examples.

Canonical inputs, read in their current state:

* ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md, Sections1--2: exact all-height
  least degree and comparison with the former rationality reduction.
* ALL_START_HEIGHT_MARKER_SAFETY.md: exact marked endpoint identity,
  original contraction normalization and all-height zero-center safety.

The retained source ledger and receipts state the exact proof inputs and\ninternal replay scope. No deposited edition is changed by this research result.

The [main checker](verify_constant_endpoint_davenport_zannier_classification.py),
[source ledger](constant_endpoint_davenport_zannier_classification_source_ledger.json)
and [receipt](constant_endpoint_davenport_zannier_classification_verification.json)
accompany this proof. A separate [Pruefer census](verify_constant_endpoint_weighted_tree_census.py)
and its [receipt](constant_endpoint_weighted_tree_census_verification.json)
check the small weighted-tree counts by labelled trees and rotation systems.
