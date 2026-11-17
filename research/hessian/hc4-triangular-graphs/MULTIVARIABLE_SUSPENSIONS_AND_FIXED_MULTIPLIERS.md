# Multivariable polynomial suspensions with arbitrary fixed multipliers

15 November 2026. The one-variable output theorem extends across
arbitrarily many polynomial variables and an arbitrary fixed polynomial
multiplier. Let `p,k>=1`, let `X=(X_1,...,X_p)` and
`Z=(Z_1,...,Z_k)` be disjoint variable sets, and take

```
A(X) in C[X], A!=0,
P(T,Z) in C[T,Z], deg_Z(P)>0,
Pi_0=(T,A(X)Y-P(T,Z)).                       (1)
```

Include `u>=0` unused affine variables and put `r=p+k+u-1`.
The ambient space has dimension `r+3`. Every polynomial closed
embedding `j:A^(r+2)->A^(r+3)` carrying a complex plane Keller pair
through (1) forces that pair to be a polynomial automorphism.
All variable counts, multiplier multiplicities, polynomial degrees,
coefficient degrees, source degrees and stable embedding degrees are
unrestricted. The multiplier is independent of `T` and of `Z`.

In particular every multivariable `xy` suspension closes. Products
of arbitrarily many independent factors, with arbitrary positive
exponents and at least one linear factor variable, are included
by taking `A` monomial. For HC4 (`r=2`), the new cases include
an arbitrary two-variable `P(T,Z_1,Z_2)` with `A=X_1`, and an
arbitrary fixed `A(X_1,X_2)` with one polynomial variable `Z`.

Claims `HC4MYIL0/1/2` are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). The conclusion concerns
plane Keller bases in every allowed stable ambient dimension. It
does not solve either plane endpoint or unrestricted HC4/JC2.

## 1. A generative polynomial transfers irreducibility (HC4MYIL0)

Put `kappa=C(T)`. For any nonconstant `P in kappa[Z]`, choose a
decomposition

```
P=R(S), S in kappa[Z], R in kappa[U],
with deg_Z(S) as small as possible.          (2)
```

Such a decomposition exists, starting with `R=U,S=P`. The positive
joint degrees multiply under composition, so minimality makes `S`
noncomposite. The explicit primary input is
[Arzhantsev--Petravchuk, arXiv:math/0608157v2, Theorem 1(i) to (iii)](https://arxiv.org/html/math/0608157v2),
now used over `kappa=C(T)`: `kappa[S]` is integrally closed in
`kappa[Z]`. The source's setup permits an arbitrary ground field.
The earlier HC4MYHX0 application used `C`; this extension of the
field scope is audited explicitly in the
[source record](sources/multivariable-suspension-import.json).

We do not need the source's special-fibre classification, Proposition 2,
or a separate generative-polynomial existence theorem. Minimal-degree
descent supplies existence, and the following field argument supplies
the precise generic geometric assertion.

### 1.1 The generic fibre of S is geometrically integral

If `h in kappa(Z)` is algebraic over `kappa(S)`, multiply it by a
nonzero polynomial in `S` clearing a monic equation's denominators.
The product is integral over `kappa[S]` and over the normal ring
`kappa[Z]`; it is therefore polynomial in `Z`, and closedness puts
it in `kappa[S]`. Thus `h in kappa(S)`. This is also the source's
Lemma 3, with the denominator step reproduced here.

Consequently `kappa(S)` is relatively algebraically closed in
`kappa(Z)`. In characteristic zero the latter field is linearly
disjoint from every finite Galois extension of `kappa(S)`: the
intersection is trivial, and restriction of the compositum Galois
group has that intersection as its fixed field, as in HC4MYHX0.
The generic coordinate ring

```
kappa(U)[Z]/(S(Z)-U)
```

embeds in `kappa(Z)`. After each finite Galois scalar extension it
still embeds in a field, and taking their union proves it integral
after an algebraic closure of `kappa(U)`. Hence `S-U` is absolutely
irreducible. This proves generic geometric integrality directly.
It makes no assertion about every special level of `S`.

### 1.2 Exact factors and source transitivity

Let `F=C(T,V)`, let `Omega` be an algebraic closure, and write
`m=deg R>=1`, with `rho` the leading coefficient of `R`. The
polynomial `R(U)+V` is separable and irreducible over `F`: it is
linear in the independent `V`, and characteristic zero excludes
generic repeated roots. Let its roots be `alpha_1,...,alpha_m`.
Then

```
P(T,Z)+V=rho product_(i=1)^m (S(T,Z)-alpha_i). (3)
```

Every `alpha_i` is transcendental over `kappa` and
`V=-R(alpha_i)`. Thus `Omega` is an algebraic closure of the
rational function field `kappa(alpha_i)`. Section 1.1 makes each
`S-alpha_i` irreducible over `Omega`; distinct roots give distinct
factors. These are exactly the geometric factors, even when `S`
is a noncoordinate polynomial with complicated fibres.

For any field extension `L/F`, the factors over an algebraic closure
of `L` are still these absolutely irreducible factors. Their Galois
permutation is exactly the permutation of the `alpha_i`. Therefore
if `R(U)+V` remains irreducible over `L`, then `P(T,Z)+V` remains
irreducible over `L`. More precisely its factor partition over `L`
is the orbit partition on those outer roots.

For a plane Keller pair `(f,g)`, the retained
[HC4MYIK3](PARAMETER_NORMS_AND_ALL_POLYNOMIAL_OUTPUTS.md#4-every-polynomial-output-in-the-plane-hc4myik3)
makes `R(f,U)+g` irreducible, since `R in C(T)[U]` is nonconstant.
It follows that

```
P(f,Z_1,...,Z_k)+g is irreducible over C(s,t)  (4)
```

for every nonconstant `P in C(T)[Z]`, every `k>=1`, and all
degrees. The transfer statement for a supplied outer irreducibility
certificate is a field statement, and also carries every retained
higher-dimensional outer certificate. The unrestricted Keller
conclusion (4) uses the plane theorem HC4MYIK3.

No full splitting-group equality is inferred. In particular the
number of geometric factors in (3) is the outer degree `m`, not
the joint polynomial degree of `P`.

## 2. The generic multiplier fibre and its class lattice (HC4MYIL1)

Factor the fixed multiplier over `C` as

```
A=c product_(j=1)^b a_j(X)^(e_j),
c in C*, e_j>=1, a_j distinct irreducible.    (5)
```

If `A` is constant, take `b=0`. All these factors are fixed under
the target Galois group. For the generic ambient fibre set

```
Q=P(T,Z)+V,
R_Omega=Omega[X,Y,Z,W]/(A(X)Y-Q(T,Z)).        (6)
```

Unused variables `W` may be absent. Then this ring is smooth and
normal, and its geometric Weil class group is

```
Cl(R_Omega)=direct_sum_(j=1)^b
                 (Z^m/Z(1,...,1)).          (7)
```

The Galois action fixes `j` and permutes the `i` coordinates in
each copy according to the outer roots `alpha_i`. The multiplicities
`e_j` do not change this lattice.

### 2.1 Smoothness and localization

First `Q` is a nonzero polynomial in variables disjoint from `X`.
The equation in (6) is primitive linear in `Y`, since `A` and
`Q` are coprime, and is therefore irreducible. To see smoothness,
work over an algebraic closure of `kappa` before adjoining `V`.
On each irreducible component of the critical locus
`partial P/partial Z_i=0`, the differential of `P` in its function
field is zero. In characteristic zero this makes `P` algebraic
over the ground field, hence constant over that algebraically
closed field. There are finitely many such components, so finitely
many critical values. None is `-V`, which is transcendental.
Thus the hypersurface `Q=0` has no critical point in its `Z`
variables after passing to `Omega`.

A singular point of (6) would have `A=0` by its `Y` derivative,
then `Q=0` by its equation, and all `Z` derivatives of `P` zero.
This is impossible. When `A` is constant the `Y` derivative is
already nonzero. Thus (6) is smooth in every case.

Inverting `A` eliminates `Y` and gives the UFD

```
R_Omega[1/A]=Omega[X,1/A,Z,W].               (8)
```

For every `j,i`, the ideal

```
D_(j,i)=(a_j(X), S(T,Z)-alpha_i)             (9)
```

is a prime divisor. The two irreducible equations concern disjoint
variable sets and are absolutely irreducible over `Omega`, so
their product coordinate ring is a domain; `Y` and `W` are free
on that divisor. These are all prime divisors removed in (8).

### 2.2 Exact relations, including repeated multiplier factors

At the generic point of `D_(j,i)`, every other `a_h` and every
other `S-alpha_l` is a unit. The equation identifies
`S-alpha_i` with a unit times `a_j^(e_j) Y`. In the normal
height-one local ring the maximal ideal is therefore generated by
`a_j`. In particular

```
div(a_j)=sum_(i=1)^m D_(j,i),               (10)
```

with multiplicity one, regardless of `e_j`.

Localizing any Weil divisor to the UFD (8) and subtracting a
principal divisor shows that the classes (9) generate. A principal
relation supported entirely there must be a unit after localization.
The only units of (8) are

```
c_0 product_j a_j^(n_j), c_0 in Omega*, n_j in Z.
```

Their divisors give exactly the independent row relations (10).
This proves (7), including torsion-freeness and the exact rank
`b(m-1)`. No class-group classification is imported.

### 2.3 The source-fixed subgroup is zero

For a transitive permutation group on `m` letters, the fixed
subgroup of `Z^m/Z(1,...,1)` is zero. Indeed if a lift `v` has
`sigma(v)-v=c_sigma(1,...,1)`, summing its entries gives
`m c_sigma=0`, hence `c_sigma=0`. Transitivity makes the lift
constant. The conclusion includes `m=1`.

The source group of a plane Keller pair acts transitively on the
outer roots by HC4MYIK3. It fixes every `j` in (7), because `A`
has constant complex coefficients. Consequently

```
Cl(R_Omega)^(Gal(Omega/C(s,t)))=0.           (11)
```

This is the property needed for the selected geometric image
component. Full target/source group equality is unnecessary.

Independence of the multiplier from `T` is an explicit hypothesis.
A parameter-dependent multiplier can have moving geometric factors
and a correlated permutation action on both indices. No conclusion
for arbitrary `A(T,X)` is inferred from (11).
For example a simultaneous swap of two multiplier factors and two
outer roots acts on the two-copy lattice as
`(v_1,v_2)->(-v_2,-v_1)` and fixes `(1,-1)`. Transitivity of the
outer action alone does not exclude this. This is a representation
control, not an asserted Keller realization.

## 3. The global carried inverse (HC4MYIL2)

Return to polynomial `P in C[T,Z]` of positive joint `Z` degree
`d`, and write `a=deg_X A`, allowing `a=0`. Let
`B_0=C[T,V]` be the global target ring. The ambient coordinate ring
has presentation

```
B_0[X,Y,Z,W]/(A(X)Y-P(T,Z)-V).               (12)
```

### 3.1 A common free filtration with no target localization

Assign the positive weights

```
weight(X_i)=d, weight(Y)=d,
weight(Z_i)=a+1, weight(W_i)=1,
weight(B_0)=0.                              (13)
```

Use a monomial order refining these weights, with `Y` first in
ties. Its relation has leading monomial `Y` times a degree-`a`
leading monomial of `A`. Its leading coefficient is a nonzero
complex constant. Monic division therefore gives the common free
`B_0`-basis of monomials not divisible by that leading monomial.
Reduction never increases (13); the bounded spans give an exhaustive
multiplicative filtration, with only `1` in degree zero.

Over the generic target field its graded equation is

```
A_a(X)Y-P_d(T,Z)=0,                         (14)
```

where both displayed leading forms are nonzero. They concern
disjoint variable sets and are coprime. Primitive linearity in `Y`
makes the associated graded ring a domain. This remains true when
`A_a` or `P_d` factors. The field decomposition (2) can have
denominators in `T`; none is used to define the global basis (12).

### 3.2 Every affine height-one target fibre is integral

Let `pi` be an irreducible in `C[T,V]`. The polynomial `P+V`
cannot become zero over `Frac(C[T,V]/(pi))`. If all positive-`Z`
coefficients vanished, `pi` would divide a nonzero polynomial
in `T`, forcing `pi=T-c`. At this prime the remaining constant
coefficient is `P(c,0)+V`, nonzero in `C(V)`. Otherwise a positive
coefficient survives.

The fixed nonzero `A(X)` stays nonzero as well, and is coprime
to this polynomial in the disjoint `Z` variables after extending
the residue field algebraically. The equation is again primitive
linear in `Y`; thus the affine fibre is geometrically integral.
This includes all leading drops and singular multiplier levels.
It does not assert integral fibres at every codimension-two target
point. Positivity of the joint `Z` degree is essential: for
`P=P(T)` the prime `V+P(T)` instead gives the reducible fibre
`A(X)Y=0` when `A` is nonconstant.

### 3.3 From zero source-fixed classes to the polynomial inverse

For a closed stable hypersurface embedding carrying a plane Keller
pair, its geometric generic image consists of disjoint affine-space
components indexed by the finite source field extension. The
component selected by the source embedding has a class fixed by
`Gal(Omega/C(s,t))`. Equation (11) makes this class zero. All its
target conjugates therefore have class zero, so every geometric
image component is principal.

At this point the retained
[HC4MYHW1 global proof](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#2-a-class-group-and-geometric-fibre-output-criterion-hc4myhw1)
applies after its principal-component step: generic degree additivity
from (14) gives the common component pencil; primitive coefficients
in the basis (12), together with every height-one integral fibre,
give global denominator saturation. Smoothness of the image and
the invertible Keller differential force its minimal label relation
to have degree one, and the retained birational Keller endpoint
gives the polynomial inverse. These are the same global steps as
[HC4MYHZ2](LEADING_FIBRE_UNIT_GATE.md#3-all-degree-nonmonic-output-equations-and-their-remaining-escape-hc4myhz2).
All their filtration, normality and affine fibre hypotheses have
been proved above for (12). No abelian class-action assumption is
needed after the source-fixed class has already been shown zero.

This proves the full assertion (1), for every `p,k>=1` and
`r>=p+k-1`, with all degrees and multiplier multiplicities unrestricted.

### Ambient and target polynomial changes

The conclusion also holds for every output pair
`Lambda composed with Pi_0 composed with Phi`, where `Phi` is
any polynomial automorphism of the ambient affine space and
`Lambda` any polynomial automorphism of the target plane. A
carried embedding becomes `Phi composed with j`, carrying the
Keller pair `Lambda^(-1) composed with K` through `Pi_0`.
The theorem makes that pair, and hence `K`, invertible. This is
an invariance of the carried-output obstruction. It does not
assert source Hessian invariance or the existence of such an
ambient normal form for an arbitrary graph embedding.

## 4. Controls and boundaries

A useful noncoordinate generative polynomial is
`S=Z_1^2+Z_2^3`. Its critical locus is zero-dimensional, whereas a
nontrivial univariate composition would contain a critical level
hypersurface. Hence it is noncomposite. Its zero fibre is singular,
so it is not a coordinate. With `R(T,U)=T U^3+U`, the generic
class lattice has rank two for each distinct multiplier factor,
although `deg_Z(P)=9`. Neither a coordinate generator nor the
joint-degree root count is legitimate here.

The exact positive controls retain actual four-variable unit Hessian.
For `A=X`, `P=(Z_1^2+Z_2^3)^3+(Z_1^2+Z_2^3)`, the carried
embedding is `(s,t,w,v)->(s,1,t+P(w,v),w,v)`. The nondegenerate
ambient quadratic `T Z_2+Y Z_1+X^2` pulls back to
`1+s v+t w+w P(w,v)`, with Hessian determinant one. Separately,
`A=1+(X_1-1)X_2`, `P=T Z^3+Z`, and the embedding
`(s,t,u,v)->(s,1,u,t+P(s,v),v)` give the same unit-Hessian
control with genuinely varying parameter coefficients. Both carried
Keller pairs are the identity; no counterexample is constructed.

The [checker](verify_multivariable_suspensions.py) and
[receipt](receipt-multivariable-suspensions.json) give exact controls
for these written proofs. The [validation](validation-multivariable-suspensions.json)
separates scoped checks from unrelated repository failures.
The [source audit](sources/multivariable-suspension-import.json)
records the exact expanded field scope of the retained closed-polynomial
input. No special-fibre classification or coordinate theorem is imported.

Arbitrary parameter-dependent or overlapping multipliers, other
ambient output geometries, obtaining compatible observables for
arbitrary HC4 graphs, actual unit-Hessian representatives in the
two endpoint source orbits, both plane endpoints and unrestricted
HC4/JC2 remain open. Formal verification, independent replay,
external review and novelty priority are not claimed.
