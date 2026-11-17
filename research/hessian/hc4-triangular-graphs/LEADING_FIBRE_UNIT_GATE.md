# Leading coefficients, localized units and a finite fibre gate

15 November 2026. Nonconstant leading coefficients need not be inverted
without control. Their exact obstruction is a localized source-unit
lattice. If the finitely many source fibres above their zeros are
irreducible, those units descend to the target parameter and the
resultant/inertia proof still closes, in every degree. Separately,
abelian root monodromy closes many nonmonic families without that
fibre assumption.

## 1. The exact localized-unit defect (HC4MYHZ0)

Let K=(f,g):A2->A2 be complex Keller, B=C[s,t], and
a(T) a nonzero polynomial. Let S be its distinct complex zeros.
The fibres f=c are smooth and reduced, because df never vanishes.
Write

\[
                  f-c=lambda_c\prod_{j=1}^{r_c} h_{cj},
                         \quad c\in S,                   \tag{1}
\]

with distinct source irreducibles. Factors from different levels are
also distinct, since f-c and f-d generate the unit ideal for c!=d.
The UFD factorization of B gives

\[
 B[1/a(f)]^*=C^*\prod_{c\in S,j} h_{cj}^{Z},\qquad
 C[f,1/a(f)]^*=C^*\prod_{c\in S}(f-c)^Z.                  \tag{2}
\]

Thus their quotient is exactly

\[
 \mathcal U_{f,S}=\bigoplus_{c\in S}
                 Z^{r_c}/Z(1,...,1),\qquad
 rank\mathcal U_{f,S}=\sum_{c\in S}(r_c-1).               \tag{3}
\]

It is torsion-free. Multiplicities of the roots of a do not affect
this group. In particular **every localized source unit descends
to a target-parameter unit if and only if all f=c, c in S, are
irreducible**. This is a global statement, not generic irreducibility.

Indeed a unit in a localization of a UFD can use only irreducible
factors of the inverted polynomial. The target units correspond to
the all-ones vector in each reduced fibre factorization (1). There
are no other relations, and each all-ones vector is primitive. This
proves (2)-(3) in unrestricted degree.

## 2. Normal monogenic covers over a vertical open target (HC4MYHZ1)

Put A=C[T,V], F=C(T,V), L=C(s,t), with T=f,V=g through K,
and A_a=A[1/a(T)]. Suppose R in A_a[Z] is monic irreducible and
Spec(A_a[Z]/R) is normal. Assume only that the finitely many
fibres f=c over c in S are irreducible. Then

\[
                             R(K,Z)\text{ is irreducible over }L.   \tag{4}
\]

All source and cover degrees are unrestricted. The assumption on
the fixed f-fibres is explicit; it is not supplied by the Keller
condition or by a silent target-coordinate change.

### Resultants descend even though they need not be constants

The etale normality argument of
[HC4MYHY0](NORMAL_MONOGENIC_COVERS.md#1-normal-finite-monogenic-covers-stay-integral-hc4myhy0)
applies on a(T)!=0. A proper monic source factorization has coefficients
in B[1/a(f)] and disjoint zero sets there, so its cross-resultant is
a unit in that ring. By (2) and the fibre assumption,

\[
               r=gamma\prod_{c\in S}(f-c)^{e_c},\qquad
                         gamma\in C^*,\ e_c\in Z.         \tag{5}
\]

It therefore belongs to F, under the target identification, and is
a unit of A_a. Every target Galois conjugate root cut has the same
r. At every target divisor outside a(T)=0, integral roots and zero
valuation of r prohibit an inertia crossing, exactly as in HC4MYHY0.
Thus those inertia groups lie in the normal core N of the proper
factor-cut stabilizer G_cut, in a finite splitting group G.

Let Q be the fixed field of N. Its finite normal target cover is
etale over U=Spec(A_a), by the same explicit ramification and purity
inputs. It need not be trivial on U. The remaining quotient is
instead a cover of the T parameter, as follows.

### Finite etale covers of the generic affine line are constant

Set k=C(T). The generic T-fibre of this quotient is finite etale
over A1_k with coordinate V. If Q/k(V) is finite Galois, then

\[
                         Q=k_0(V),\qquad k_0/k
                         \text{finite Galois}.             \tag{6}
\]

Here is the needed descent proof. Let k0 be the algebraic closure
of k in Q. It is finite: degrees of its finite subextensions are
bounded by [Q:k(V)], so an increasing chain stabilizes. It is Galois
over k because Q/k(V) is normal. The normalization of A1 over k0
in Q remains finite etale, by constant-field base change and selection
of the field component. Take the smooth projective curve with function
field Q over k0. It is geometrically integral in characteristic zero,
since k0 is relatively algebraically closed in Q. After extending
constants to an algebraic closure, it gives a connected degree-e
cover of P1, unramified away from infinity.

By the explicitly imported characteristic-zero
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B),

\[
                  2g-2=-2e+\sum_{p\mid infinity}(e_p-1)
                       =-e-r,                            \tag{7}
\]

where r>=1 is the number of points above infinity. Since g>=0,
e+r<=2; hence e=r=1, proving (6). The finite normalization of
the affine V-line includes every point above finite V; this is why
there are no omitted finite boundary points in (7).

### The regular source parameter eliminates that quotient

The earlier [HC4MYHX0](REGULAR_PARAMETER_MONODROMY.md#1-nowhere-critical-polynomial-parameters-are-regular-hc4myhx0)
makes C(f) relatively algebraically closed in L and preserves every
finite Galois group over that parameter. Thus the source stabilizer
H=Gal(Omega/L) maps onto Gal(k0/k), which by (6) is the entire
quotient G/N. But H's image in G preserves the selected source
factor, so lies in G_cut, and N lies in G_cut as well. Surjectivity
then gives G=H*N subset G_cut, impossible for a proper cut of an
irreducible polynomial. This proves (4). It asserts root irreducibility,
not equality of full splitting groups after source base change.

There is also a finite-normalization consequence. If F subset E
subset L is finite and its normalization over A_a is globally
monogenic there, A_{a,E}=A_a[h], then E=F under the same finite
fibre hypothesis. The integral h belongs to B[1/a(f)], and its
minimal monic normal equation would be irreducible over L by (4)
while having the root h in L. Thus it has degree one. A nontrivial
intermediate normalization becoming monogenic after only a vertical
target localization therefore forces at least one reducible f-fibre
among the inverted levels. A generic primitive element still does
not supply a normal monogenic presentation.

## 3. All-degree nonmonic output equations and their remaining escape (HC4MYHZ2)

Let P(T,z) be any polynomial of positive z degree m, with arbitrary
nonzero leading coefficient a(T). For r>=1, consider a genuine closed
embedding j:A^(r+2)->A^(r+3) carrying a plane Keller pair K through

\[
                            Pi=(T,xy-P(T,z)).              \tag{8}
\]

Then **K is invertible if every K1=c over a zero of a is irreducible**.
Only those finitely many full source-base fibres are required. All
T, z, leading-coefficient and embedding degrees are unrestricted.

Over a(T)!=0, the root equation (P(T,Z)+V)/a(T) is monic and
its total cover is Spec C[T,1/a(T),Z], a smooth normal surface.
Apply (4) to obtain a transitive source action on all its generic
roots. The root-divisor class lattice Z^m/Z(1,...,1) has zero
invariants for every transitive permutation subgroup, as proved in
HC4MYHY1. Thus every source-component divisor class is zero.

The global denominator step still works across leading-coefficient
zeros; no localization is left in the embedding conclusion. Use
weights m,m,2 for x,y,z and degree zero on A. The monic xy relation
gives the common free A-basis. Its generic graded equation
xy=a(T)*z^m is a domain because a(T) is nonzero generically.

At any affine height-one target prime pi, P(T,z)+V cannot become
the zero polynomial. If all positive-z coefficients vanished modulo
pi, pi would divide a nonzero polynomial in T, so pi=T-c up to
scalar; then the constant coefficient is P(c,0)+V, nonzero in C(V).
Otherwise a positive-z coefficient survives. In every case
xy-(P(T,z)+V) is geometrically irreducible, by the primitive linear
polynomial in y with coefficients x and a nonzero polynomial in z.
All affine height-one fibres are therefore geometrically integral,
even when their graded rings degenerate or their z degree drops.
The geometric generic ring is smooth, since V is transcendental over
C(T) and P(T,z)+V has no repeated generic root in characteristic zero.
Localizing x gives the same root-divisor class calculation as before;
the nonzero generic leading coefficient is a field unit. Integral
height-one fibres do not assert integrality of every closed target
fibre: at a degree-drop point the latter can still have xy=0.
The component-pencil, global primitive saturation and smooth minimal
label argument of HC4MYHW1 give the polynomial inverse for K.

The same global filtration proof gives a second, independent closure:
**if the generic root group of P(T,z)+V over C(T,V) is abelian,
no irreducibility assumption on the special K1-fibres is needed**.
Its generic root-divisor class group is torsion-free with abelian
action, so HC4MYHW1 applies directly. Consequently:

- Every P of z degree one or two is closed, with arbitrary T coefficients
  and arbitrary nonconstant leading coefficient.
- Every P=a(T)*(z-b(T))^m+c(T), with a nonzero and a,b,c in C[T],
  is closed for every m>=1 without a special-fibre assumption. The
  roots split in a cyclic degree-m extension of C(T,V): after translating
  z the m-th power equals -(V+c(T))/a(T), whose valuation at V+c(T)
  is one. Eisenstein proves its degree and C supplies all roots of unity.

Thus a hypothetical nonautomorphic realization of (8) must have
**both nonabelian generic root monodromy and a reducible full K1-fibre
at a zero of the leading coefficient**. Equivalently it must have
positive localized-unit defect (3) there. Nonconstant a by itself,
a root pole by itself, or increasing z degree is not an escape.

For genuine positive controls, j(s,t,w)=(s,1,t+P(s,w),w) is closed
with inverse (T,y-P(T,z),z), for every P here, including leading-degree
drops. Add an unused v. The nondegenerate quadratic T*v+y*z+x^2
has actual pullback 1+s*v+t*w+w*P(s,w), with unit Hessian and the
same explicit gradient inverse from HC4MYHY1. These controls have
invertible gradient; no nonautomorphic Keller pair is constructed.

## 4. Why the fibre condition is real

The polynomial f=s+s^2*t has no critical points, by
(1-2st)f_s+4t^2*f_t=1, but f=0 consists of the disjoint components
s=0 and 1+st=0. On B[1/f], both s and 1+st are units. Their
classes modulo C[f,1/f]^* give the defect Z^2/Z(1,1), of rank one.
The unit s does not belong to C(f): indeed L=C(f,s) with
t=(f-s)/s^2 and s transcendental over C(f). Thus generic parameter
regularity does not give localized unit descent when special fibres
are reducible. This is a polynomial submersion control, not a Keller
counterexample or a stipulated Keller component with a mate.

No theorem that all fibres of the fixed Keller component K1 are
irreducible is used. Nor is a target-coordinate reduction silently
applied to the fixed ambient T. Such a change would also change
the two ambient observables and needs its own verified normal-form
certificate. The remaining positive-defect, nonabelian branch, both
plane endpoints and unrestricted HC4/JC2 remain open.

## 5. A source-audit boundary: affine etale is not finite

The [source record](sources/leading-fibre-import.json) retains
Tian--Sun, AIMS Mathematics 9 (2024), 30406-30412, but **imports
none of its claimed Keller conclusions**. In the proof of its
Theorem 2.8 the boundary of a completed affine source curve is
identified with the fibre over infinity from etaleness on the affine
part. This is not justified without excluding finite boundary values.

The exact control Gm->A1, u->u, is etale. On its completion P1,
the omitted boundary is {0,infinity}, whereas the fibre over infinity
is only {infinity}. The same issue can occur after resolving a
polynomial fibre at infinity: omitted points may map to finite values
when the affine map is not finite. Our use of (7) instead starts
with an actual finite etale cover of the whole affine line, so that
all points over finite V are included by its finite normalization.
This audit does not use the disputed shortcut to resolve JC2.

The [checker](verify_leading_fibre.py) and
[receipt](receipt-leading-fibre.json) give exact finite controls,
separate from these written all-degree proofs and the explicit
Riemann--Hurwitz input. No formal verification, independent replay,
external review or novelty priority is asserted.
