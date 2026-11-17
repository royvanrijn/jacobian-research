# A rank-two pencil lattice removes the global-frame hypothesis

15 November 2026. The global denominator step needs a primitive
pencil coordinate, rather than a common free basis of the entire
ambient ring. Faithful flatness and geometrically integral height-one
fibres produce that coordinate over the plane target. Consequently
every disjoint moving multiplier suspension satisfying the exact
content-fibre condition closes, with no leading-coefficient or
global-frame certificate. All degrees and admitted stable dimensions
are unrestricted. The unrestricted HC4/JC2 problem and both plane
endpoints remain open.

Claims HC4MYIN0/1 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). Previous proof notes and
receipts are preserved; their sufficient global frames remain valid.

## 1. A primitive global coordinate for any generic pencil

Let `A=C[T,V]`, `F=Frac(A)`, and let `R` be a normal finitely
generated complex domain, faithfully flat as an `A`-algebra. Suppose
each height-one generic fibre

```
R tensor_A kappa(pi),  kappa(pi)=Frac(A/(pi)),
```

is geometrically integral. Let `h0 in R_F=R tensor_A F` be
transcendental over `F`. Then there exists `h in R` with

```
F+F h = F+F h0,
(F+F h0) intersection R = A direct_sum A h,
F[h] intersection R = A[h].                         (1)
```

No filtration of `R` over `A` is assumed. A generic filtration
will be used only later, to construct the pencil from the image
components. This lemma descends any pencil already supplied.

### 1.1 Pullbacks of target prime divisors are prime

Faithful flatness gives `IR intersection A=I` for every ideal
`I` of `A`, and hence `R intersection F=A`. For an irreducible
`pi in A`, the ring `R/pi R` is flat, hence torsion-free, over
the domain `A/(pi)`. It embeds in its generic fibre, which is
a domain by hypothesis. Thus `pi R` is a proper prime ideal.
It has height one, and `R_(pi R)` is a DVR by normality.
Its uniformizer is `pi` itself. In particular, for every `q in F*`,

```
valuation_(pi R)(q)=valuation_pi(q).                (2)
```

Distinct target primes give distinct primes of `R` by contraction.
This argument excludes hidden components of a pulled-back divisor;
generic integrality without flatness would not suffice.

### 1.2 One derivation bounds every coefficient denominator

Multiply `h0` by a nonzero element of `A` to arrange `h0 in R`;
its two-dimensional `F`-span is unchanged. Put

```
W=F+F h0,  N=W intersection R.                     (3)
```

There is an `A`-derivation `D:R->R` with `d=D(h0)!=0`.
Indeed extend `h0` to a transcendence basis of `Frac(R)/F`,
differentiate with respect to `h0`, and extend through the finite
separable algebraic extension in characteristic zero. Multiplication
by a common denominator of the derivatives of finitely many algebra
generators gives the required derivation of `R`.

If `q+p h0 in N`, with `q,p in F`, then

```
D(q+p h0)=p d in R.                                (4)
```

Only finitely many prime divisors of `R` occur in `div(d)`.
Let `a in A` be the product of the target primes whose pulled-back
prime divides `d`, each with exponent `valuation_(pi R)(d)`.
Equation (2) applied to (4) gives `p in a^(-1) A`. Then
`q in a^(-1) R intersection F=a^(-1) A`. Consequently

```
A+A h0 subset N subset a^(-1)(A+A h0).              (5)
```

Since `A` is Noetherian, `N` is a finite torsion-free module of
rank two. This is an actual uniform denominator bound, not an
assumption that a finite-type flat algebra is free as a module.

### 1.3 Reflexivity and the nowhere-vanishing constant section

Inside `W`,

```
N = intersection_(height(pi)=1) N_pi.               (6)
```

To verify the reverse inclusion, write an element as `r/b` with
`r in R,b in A`. At every height-one prime of `R` above a factor
of `b`, membership in `N_pi` removes the possible pole. The
principal prime calculation in Section 1.1 shows these are all
prime divisors of `b`. At every other prime divisor the denominator
is a unit. Normality of `R` now puts `r/b` in `R`.

For clarity, the standard double-dual criterion in this instance
follows directly. For a finite torsion-free `A`-module, its
double dual is the intersection of its height-one localizations:
one direction uses that finite torsion-free modules over a DVR
are free; in the other direction every linear functional has value
in every `A_pi`, hence in `A`. Thus (6) says `N` is reflexive.

The explicitly imported [Stacks 0B3N](https://stacks.math.columbia.edu/tag/0B3N)
states that a coherent reflexive module on a regular scheme of
dimension at most two is finite locally free. It makes `N` locally
free of rank two on `Spec A`. We do not claim this for regular
bases of arbitrary dimension.

The element `1 in N` is nonzero in every closed residue fibre of
`N`: otherwise `1 in m N subset m R` for a maximal ideal `m`
of `A`, contradicting faithful flatness. Thus it is a locally
unimodular section, and `N/A*1` is an invertible module. Every
invertible module over the UFD `C[T,V]` is free: view it as an
invertible fractional ideal and remove its finite height-one divisor;
the resulting reflexive rank-one ideal is `A`. Therefore

```
N=A*1 direct_sum A*h                               (7)
```

for some `h in R`. This uses neither a global basis for `R` nor
Quillen--Suslin for a higher-rank module.

### 1.4 Primitive reductions and polynomial saturation

For every height-one `pi`, the reduction of `h` in
`R tensor_A kappa(pi)` is not a scalar. Otherwise, for some
`c in A_pi`, the element `(h-c)/pi` would belong to
`W intersection R_pi=N_pi`, contradicting the basis `(1,h)` in
(7). Geometric integrality in characteristic zero makes every
nonscalar element transcendental over `kappa(pi)`: an algebraic
separable minimal polynomial would split into linear factors in
a geometrically integral scalar extension, forcing the element
to be a scalar already by linear independence.

If `B in F[S]` has negative minimum coefficient valuation at
`pi`, multiply it by the compensating power of `pi`. Its reduction
is a nonzero polynomial and cannot vanish at this transcendental
reduction of `h`. Thus `B(h)` cannot be in `R_pi`. Applying this
at every height-one target prime proves `F[h] intersection R=A[h]`,
and finishes (1).

## 2. An output criterion requiring only a generic filtration (HC4MYIN0)

Let `K:A2->A2` be plane Keller, and let a genuine polynomial
closed embedding `j:A^(r+2)->A^(r+3)`, `r>=0`, carry `K` through
ambient observables `Pi=(P1,P2)`. Set `A=C[P1,P2]` and
`R=C[ambient variables]`, identifying `A` with `C[T,V]`.
Assume:

1. `R` is faithfully flat over `A`, and every height-one generic
   fibre is geometrically integral.
2. `R_F` has an exhaustive multiplicative nonnegative filtration
   with finite-dimensional bounded pieces and degree-zero part `F`.
   After extension to `Omega=F_bar`, its associated graded ring
   is a domain with degree-zero part `Omega`.
3. `R_Omega` is normal and every geometric generic image component
   of `j` has class zero. It suffices that
   `Cl(R_Omega)^Gal(Omega/L)=0`, where `L=C(s,t)` is the selected
   plane source field in `Omega`.

Then `K` is a polynomial automorphism, in every degree. Here a
generic filtration is sufficient; no common global free filtration
or global monic relation is required.

The prime image equation `E` has geometric generic zero set a
disjoint union of `[L:F]` copies of `A^r`. Condition 3 makes these
components principal. The generic degree-additivity and normalized
component-pencil argument of
[HC4MYHV0](SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md#1-a-filtered-output-base-criterion-hc4myhv0)
uses only an `F`-basis compatible with condition 2: conjugate
component generators have equal degree, disjointness makes each
a nonzero constant modulo every other, and equal degrees make
their quotient scalar. Normalizing the last basis coefficient
gives generators `h0-alpha_i`, with `h0 in R_F` by Galois descent.
Positive generic degree makes `h0` transcendental over `F`.
The residual generic unit in `E` is a scalar by condition 2.

Section 1 supplies `h in R` with `F[h]=F[h0]` and the literal
global identity

```
E=B(T,V,h), B in A[S], deg_S B=[L:F].               (8)
```

There is no denominator left in (8). As in HC4MYHV0, `B` is
irreducible: substitution is injective and every nonconstant factor
substitutes to a nonconstant element of the polynomial ring `R`.
On the image, `h composed with j=a(s,t)` is independent of the
stable variables, since it is a constant label on each geometric
generic image component. This gives its minimal relation `B`.

Smoothness of the image and invertibility of `DK` imply
`B_S(K,a)` never vanishes. It is therefore a nonzero constant `c`.
The polynomial `B_S-c`, of smaller label degree, vanishes at the
same label. Minimality forces it identically zero, so `B=c S+b`
and `[L:F]=1`. The retained birational Keller endpoint gives the
polynomial inverse. This argument imposes no Hessian condition
and does not assert that arbitrary observables satisfy conditions 1--3.

## 3. All exact-content disjoint moving suspensions (HC4MYIN1)

Take arbitrary `A(T,X) in C[T,X]`, `A!=0`, and
`P(T,Z) in C[T,Z]` of positive joint `Z` degree. The sets `X,Z`
are disjoint, with `p=#X>=1,k=#Z>=1`. Let `delta_A(T)` be the
gcd of all `X` coefficients of `A`. Suppose that for every zero
`c` of `delta_A`,

```
P(c,Z) is nonconstant and closed over C.            (9)
```

Then every genuine closed stable embedding carrying a plane Keller
pair through

```
Pi=(T,A(T,X)Y-P(T,Z)), r>=p+k-1,                    (10)
```

forces that pair invertible. All parameter, coefficient, source,
embedding and variable degrees and all admitted stable dimensions
are unrestricted. No global-frame, leading-unit, leading-divisibility,
Bezout-template or special source-fibre hypothesis remains.
In particular **every primitive multiplier** (`delta_A` constant)
is covered, with arbitrary `P`. Condition (9) may also allow
nonconstant multiplier content.

Here is the verification of every hypothesis of Section 2.

- **Faithful flatness.** The ambient coordinate ring is the polynomial
  ring in `T,X,Y,Z` and unused variables. At `T=c` its second
  observable is nonconstant: either `A(c,X)!=0`, or (9) makes
  `P(c,Z)` nonconstant. Hence `T-c,V-d` is a regular sequence in
  every local ambient ring above `(c,d)`, and every closed target
  fibre is nonempty of the expected dimension. At target primes of
  height one, torsion-freeness over the DVR already gives flatness;
  at height zero it is automatic. At height two the explicitly
  imported [miracle flatness lemma, Stacks 00R4](https://stacks.math.columbia.edu/tag/00R4)
  applies to the regular target and regular ambient local rings,
  with dimension drop exactly two. Thus the morphism is flat.
  Every maximal target ideal has proper extension, so it is faithfully
  flat. This proves surjectivity of this ambient observable map,
  not surjectivity of a hypothetical carried Keller map.
- **Height-one fibres.** The exact criterion HC4MYIM2 says precisely
  that all these affine generic fibres are geometrically integral
  if and only if (9) holds. It uses disjointness and primitive
  linearity in `Y` away from the content levels, and the generative
  polynomial theorem at those levels.
- **Generic filtration.** Write `a=deg_X A>=0,d=deg_Z P>0`.
  Over `F` use weights `d` on `X,Y`, `a+1` on `Z`, and one
  on unused variables. Over a field every nonzero leading coefficient
  is a unit, so standard monomials give the filtered `F`-basis.
  The geometric generic graded equation is
  `A_a(T,X)Y-P_d(T,Z)`. It is primitive linear in `Y` over
  `Omega`, since its coefficients have disjoint variables, and
  is a domain. These are only generic assertions; no global basis
  is being inferred from them.
- **Principal image components.** HC4MYIM1 gives geometric generic
  smoothness, the exact moving multiplier class lattice, and zero
  source-fixed class subgroup for every such disjoint `A,P`.
  That theorem uses HC4MYIM0 parameter-stable outer irreducibility.

Thus HC4MYIN0 applies. Ambient and target polynomial automorphism
orbits inherit the conclusion, exactly as in HC4MYIL2. No nonlinear
source Hessian invariance is used.

For one `Z`, (9) means that every nonconstant content level has
degree one. With multiple `Z` variables it permits any nonconstant
closed polynomial, including products such as `Z1*Z2`, whose zero
level can be reducible although its generic level is integral.
The uniform example `A=T(1+(T+1)X)`, `P=T Z^m+Z`, every
`m>=2`, satisfies (9), including the content zero at `T=0`.
The theorem does not stipulate that any specified nonautomorphic
pair can be embedded through this example.

A genuine HC4 positive control is available for arbitrary
`P(T,Z)` and `A=1+T X1^2+(T+1)X1 X2^2`:

```
j(s,t,u,v)=(s,0,u,t+P(s,v),v),
Pi composed with j=(s,t),
Q=T X2+Y Z+X1^2,
Q composed with j=s u+t v+v P(s,v).
```

The inverse on the image reads `(T,Y-P(T,Z),X2,Z)`, and the actual
four-variable Hessian determinant of the pullback is one. This
checks compatibility with actual closed embeddings and a unit-Hessian
potential; its gradient is invertible, and it constructs no
nonautomorphic Keller pair.

## 4. Why the hypotheses still matter

Faithful flatness cannot be replaced merely by a flat map with
integral height-one generic fibres. Take

```
R=A[x,y]/(T x+V y-1).
```

It is a smooth normal domain, flat over `A`, with geometrically
integral height-one generic fibres. Its image omits `(0,0)`, so
it is not faithfully flat. Generically it is an affine line.
For `W=F+F x`, the lattice `W intersection R` is
`A x direct_sum A y`, with constant section `1=T x+V y`.
This section vanishes as a vector at the omitted base point;
the quotient is the noninvertible ideal `(T,V)`. Thus it has no
basis `(1,h)` as in (7). The omitted point is exactly where the
Section 1.3 argument would fail. This is an algebraic control,
not a polynomial-ambient Keller counterexample.

Geometric height-one integrality remains separate: the retained
`(T,T Y+Z^2)` denominator trap is faithfully flat but its
`T=0` generic fibre splits. It still has
`(Z^2-V)/T=-Y`, and does not satisfy this theorem.

If `A(c,X)=0` and `P(c,Z)` is constant, HC4MYIM2 already forbids
any carried Keller pair by Jacobian divisibility. Nonconstant
composite content levels remain outside the integral-fibre theorem.
Failure there is not asserted to produce a Keller realization.
Overlapping multiplier/root variables, other output geometries,
compatible observables for arbitrary HC4 graphs, actual endpoint
representatives, both plane endpoints and unrestricted HC4/JC2
remain open. A missing global frame is no longer a residual in
the disjoint exact-content family.

The [source audit](sources/reflexive-pencil-import.json) explicitly
imports only the reflexive-surface and miracle-flatness lemmas.
The [checker](verify_reflexive_pencil.py),
[receipt](receipt-reflexive-pencil.json) and
[scoped validation](validation-reflexive-pencil.json) distinguish
exact controls and byte pins from the written general proof.
Formal verification, independent replay, external review and novelty
priority are not claimed.
