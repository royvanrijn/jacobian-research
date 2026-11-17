# Canonical boundary and stable complete-intersection obstructions

Claims HC4MYKN0/1/2. [MATH_STATUS.json](../../MATH_STATUS.json) is the
authority. This returns to the global finite normalization of a Keller
extension. It excludes torsion-canonical intermediate normalizations and
every global complete-intersection presentation,with arbitrarily many
helpers and polynomial or Laurent stabilizations. The criterion does NOT
prove such a presentation always exists. Full HC4/JC2 remains open.

The full-normalization class-group calculation is already retained in
[KRP1,Section5](../../foundations/verified/KELLER_RITT_PRODUCT_THEOREM.md).
The new use here is its canonical ramification class and an argument
for EVERY embedded intermediate field,where no intermediate affine-space
model is assumed. The degree-six smooth-cover discussion of
[finite-flat trace lattices,Section6](https://github.com/alok/jacobian-two/blob/main/docs/finite-flat-trace-lattices.md)
already identifies the analogous non-Frobenius and global-presentation
obstructions in that restricted setting. We credit that source;we do
not claim the degree-six observation or the localization sequence anew.

## 1. Global normalization and its boundary: HC4MYKN0

Let K:X=A^n_C -> Y=A^n_C be polynomial Keller,n>=1. Set

    B=C[T1,...,Tn], R=C[z1,...,zn],
    F=Frac(B), L=Frac(R), B -> R through K.

The finite normalization S of B in L is a subring of R:an element
integral over B is integral over the normal ring R and hence lies in R.
Its spectrum Z is finite over Y,and the canonical lift X->Z is an
open immersion by the retained Zariski--Main argument. It is birational
and quasi-finite because K is. This is the FULL source field L;
no such birational assertion is made for a proper intermediate field.

Write D1,...,Ds for the prime divisors of Z outside X. The retained
KRP1 localization argument gives

    Cl(Z)=direct_sum_(i=1)^s Z[D_i].                           (1)

For completeness,the proof is short. Any Weil divisor restricts to a
principal divisor on the factorial X;subtract that principal divisor
to leave only the D_i. A relation sum a_i D_i=div_Z(f) makes f a
rational function with no zero or pole on X. Normality makes it an
actual unit on X,and R*=C*,so f is constant and all a_i vanish.
This proves both surjectivity and independence. Possible omitted
higher-codimension points do not alter this calculation.

Let pi:Z->Y be the finite map. At a source prime D over a target prime
divisor,write e_D for its ramification index. Since the characteristic
is zero,the pulled-back target volume form gives the canonical divisor

    R_pi=sum_D (e_D-1)D,
    [K_Z]=[R_pi].                                             (2)

This is a Weil canonical divisor on the NORMAL Z;Z need not be smooth,
Gorenstein or globally a complete intersection. Only its height-one
local rings are used. Here is the local differential calculation. At
generic source and target prime points,choose target transverse parameter
t and separating tangential parameters y2,...,yn. In the source DVR,
t=a*u^e with a a unit. The tangential differentials stay independent
because the residue extension is separable. The wedge of

    dt=e*a*u^(e-1)du+u^e da

with those tangential differentials has order exactly e-1:its leading
coefficient is a unit and the second summand has strictly higher order.
This is the usual tame codimension-one calculation;the inspected
[Stacks0C1B](https://stacks.math.columbia.edu/tag/0C1B) curve calculation
uses the same local formula. No wild-characteristic formula is inferred.

K is etale on X,so every ramified D is among the D_i. Therefore (1)--(2)
give the exact free-basis coordinates

    [K_Z]=sum_(i=1)^s (e_i-1)[D_i].                            (3)

If [L:F]>1,at least one coefficient is positive. Otherwise pi is
unramified in codimension one;retained purity makes it etale everywhere,
and the retained triviality of connected finite etale covers of complex
affine space gives degree one. Thus a hypothetical nonautomorphic plane
Keller map has a NONZERO,NONTORSION canonical class in a free boundary
class group. This is a necessary global obstruction,not a contradiction:
the class is allowed to be nonzero.

The source open cannot be recovered by inverting one element of S,or
any finite list of elements. Such a list reduces to its product f. If
X=D(f) inside Z,then f restricts to a unit on X and hence is constant
by dominance. It cannot remove a nonempty boundary. This is another
description of why localizing one primitive equation does not restore
the entire affine-space source.

## 2. Every embedded intermediate normalization: HC4MYKN1

Now take ANY embedded finite intermediate field

    F subset E subset L,

and let S_E be the finite normalization of B in E. It still lies in
R,and gives a dominant quasi-finite map X->Z_E=Spec S_E followed by
the finite pi_E:Z_E->Y. Do NOT assume X->Z_E is birational or Z_E
has an affine-space reconstruction open. Those are additional properties
of selected polynomial factorizations in KRP1,not properties of an
arbitrary intermediate field.

We claim

    [K_(Z_E)] torsion in Cl(Z_E)  ==>  E=F.                    (4)

Use the same effective ramification divisor R_E representing K_(Z_E).
If its class is torsion,there exist m>=1 and a nonzero rational function
g on Z_E with

    div_(Z_E)(g)=m R_E.                                       (5)

There are no codimension-one poles,so normal affine Hartogs makes
g an element of S_E. It is consequently a source polynomial in R.
If this polynomial is nonconstant,it has an irreducible source zero
divisor C. Quasi-finiteness of K makes C map onto a prime divisor D
of Z_E over a target prime P. By (5),D occurs in R_E,so e(D/P)>1.
But multiplication of ramification indices in the actual embedded
valuation tower gives

    e(C/P)=e(C/D)*e(D/P),                                    (6)

where the left side is one because K is etale. This is impossible.
Thus g is a source unit,so it is a nonzero complex constant. Injectivity
S_E->R then makes g that same constant in S_E. Equation (5) gives
R_E=0. Purity and affine simple connectedness again make pi_E a trivial
connected finite etale cover. Hence E=F,proving (4).

This proof uses effectiveness,scalar source units and the ACTUAL
intermediate valuation tower. It requires neither a global monogenic
presentation nor an affine-space model of E. It also works for an
etale dominant generically finite morphism from a normal affine source
with only scalar units;factoriality was used in (1),not in (4).

For E=L,(4) gives an inverse criterion for K. Once L=F,the retained
Keller denominator-removal argument makes the rational source inverse
polynomial. For a proper E,(4) only excludes that intermediate field;
it does not assert the whole extension L/F trivial.

As a consistency check,a finite Galois normalization over affine space
has torsion canonical class:ramification indices above a given target
prime P are all the same e_P,and

    R_E=sum_P ((e_P-1)/e_P)*pi_E^*(P).

Clearing denominators makes it principal,as every target prime is
principal. Thus (4) reproduces the already-retained KMFI1 exclusion of
nontrivial Galois intermediate fields. That Galois exclusion is not a
new result here. The additional target is nongalois normalizations
whose canonical class happens to be torsion.

## 3. Arbitrarily many helpers and stable presentations: HC4MYKN2

Suppose S_E has ANY global complete-intersection presentation as an
affine complex algebra. A Koszul resolution gives

    omega_(S_E)=Ext^c_P(S_E,omega_P) isomorphic to S_E,

where P is the ambient polynomial ring and c is the regular-sequence
length. The last Koszul module and ambient top differential module
both have global bases;the top dual cohomology is exactly S_E. Thus
the canonical class is zero,and (4) forces E=F.

In particular NO nontrivial intermediate normalization can have a
global square presentation

    S_E = B[U1,...,Ur]/(f1,...,fr),
    (f1,...,fr) a regular sequence,                            (7)

for ANY helper number r or ANY polynomial degrees. This strictly
extends the earlier retained HY2 one-integral-generator obstruction
as a presentation criterion;it does not prove every multi-helper
normalization has form (7).

Adding trivial variables cannot rescue a global complete intersection.
If

    S_E[t1,...,ta,s1^+/-1,...,sb^+/-1]

had a global complete-intersection presentation,quotient successively
by the regular sequence t1,...,ta,s1-1,...,sb-1. Lifting these elements
to the ambient polynomial presentation appends them to its regular
sequence,and gives a global complete-intersection presentation of S_E
itself. This contradicts (4). Hence (7) is excluded after EVERY finite
polynomial or Laurent stabilization as well. Nontrivial vector-bundle
total spaces or nonstandard open modifications are not these products.

The canonical module of S_E is a finite torsion-free module of rank
one. If it had one GLOBAL generator,the surjection S_E->omega_(S_E)
would have zero kernel by generic rank one and the domain property.
It would be free,and (4) would apply. Thus a nontrivial intermediate
normalization needs at least TWO global canonical-module generators.
This does NOT say that its local Cohen--Macaulay type is at least two.
It can be smooth or locally Gorenstein,with an invertible but globally
nontrivial canonical line.

For a finite normal surface algebra over B=C[T1,T2],finite duality
identifies its canonical module with Hom_B(S_E,B). The inspected
[Stacks0FKW](https://stacks.math.columbia.edu/tag/0FKW) finite dualizing-
module formula is the one new external proved input used for this
identification. Equivalently its trace-dual fractional lattice is the
inverse different:at a tame prime of index e it has order -(e-1),
so its associated Weil class is the POSITIVE ramification class in
(2),not its negative. Consequently S_E cannot be a globally Frobenius
B-algebra,which would require Hom_B(S_E,B) isomorphic to S_E.
Freeness over B is not freeness as an S_E-line.

One can check the sign and the canonical identification directly in
codimension one. Over a tame completed DVR model t=u^e,the trace dual
of the basis 1,u,...,u^(e-1) is the fractional module u^(1-e)S/e.
Unramified residue-field extension contributes a unit. Intersecting
these lattices over all height-one primes is valid because B and S_E
are normal:the trace condition is tested at all height-one primes of
B,and the resulting S_E-module is reflexive. It is exactly the module
of rational top forms for the divisor R_E in (2). Thus the relative
finite-duality formula is being matched to the canonical class,not
mistaken for a free module just because the base canonical bundle is
trivial.

If the full plane normalization is presented using two extra algebra
generators,its defining ideal has height two. Two defining equations
would form a regular sequence in the regular polynomial ambient ring,
giving (7). Under nonautomorphy it therefore needs at least THREE
defining relations in every such two-helper presentation. With any
number of helpers the obstruction remains global canonical triviality.

### A genuine non-Keller control

The smooth affine quadric

    Z: X*Y=Z0^2-1

is a global complete intersection with trivial canonical class. It has
the actual affine-plane open chart

    (s,t) |-> (X=s,Y=2t+s*t^2,Z0=1+s*t).                       (8)

Its complement is the line X=0,Z0=-1. The inverse t=(Z0-1)/X on
D(X) and t=Y/(Z0+1) on D(Z0+1) glues;these opens cover precisely
that chart. Projection (X,Y) is finite of degree two. The composite
map on the affine-plane chart is

    (s,t) |-> (s,2t+s*t^2), J=2(1+s*t).                       (9)

It is not Keller on the whole chart. The ramification Z0=0 meets
the chart in the actual hyperbola 1+s*t=0. Thus a normal finite
complete-intersection cover containing an affine-plane open is not
alone a contradiction. The etale source hypothesis that places ALL
ramification outside that open is indispensable.

Characteristic zero is indispensable too. In characteristic p,the map
(x,y)->(x^p-x,y) is finite etale of degree p;its source normalization
is an affine plane,globally monogenic and canonically trivial. Finite
etale covers of affine space need not be trivial there. This is a
positive-characteristic control,not an excluded example of (4).

## 4. Relation to the actual remaining HC4 gaps

For ANY actual four-variable unit-Hessian potential h,apply the theorem
directly to the Keller map grad(h) in dimension four. If its finite target
normalization has torsion Weil canonical class,the gradient has a polynomial
inverse. A hypothetical failure therefore has a nontorsion canonical
normalization and no global complete-intersection presentation,even stably.
This direct HC4 criterion uses NO Tang classification or native graph
alignment. The plane specialization is useful through the retained HC4/JC2
reduction,but is not needed for this implication.

The primitive-helper image in IV1/IV2 is usually a NONNORMAL order.
It is not the finite NORMALIZATION in (4). A primitive field generator
or a globally embedded HC4 hypersurface must not be identified with
S_E without proving equality and finiteness over the Keller outputs.
Nor does a local complete-intersection chart or smoothness trivialize
the canonical module globally. The earlier degree-six trace-lattice
source explicitly preserves this distinction;its proposed surviving
classes are nontorsion and are compatible with the present obstruction.
The archived quartic QCL4/QLD4 packets likewise retain a nontrivial
canonical/different class and are not excluded by merely restating it.

What is now excluded is an entire all-degree construction architecture:
normal embedded finite covers with torsion canonical class,including
global complete intersections and their trivial stabilizations. A useful
next result must either prove the full normalization has such a property
from additional actual global Keller data,or control the remaining
nontorsion canonical module and multiplication algebra. No such automatic
trivialization is supplied. A proper intermediate field with nontorsion
canonical class also need not admit any polynomial intermediate model.
Full HC4/JC2,the non-CI plane normalization case,JA/JB global rational
recovery,native critical non-CI recognition and the endpoint inversions
remain open.

The [checker](verify_canonical_boundary.py),
[receipt](receipt-canonical-boundary.json) and
[context](sources/canonical-boundary-context.json) bind valuation,
canonical-class,regular-sequence and chart controls with unchanged
retained evidence. The arbitrary-field-extension argument and divisor
proofs are written mathematics;finite checks do not certify purity or
the geometric classification. No old registered mathematical checker
was replayed. Independent replay,formal verification and external review
are not claimed.
