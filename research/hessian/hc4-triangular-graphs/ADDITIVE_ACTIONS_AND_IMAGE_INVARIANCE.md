# Additive actions replace positive filtrations or force image invariance

15 November 2026. A locally nilpotent derivation moving the image
equation supplies the generic component pencil. Its degree-zero
ring may be large; neither finite-dimensional filtered pieces nor
a positive degree on every nonconstant function is needed.
Consequently a hypothetical nonautomorphic carried plane Keller
pair must have its image equation fixed by every relative additive
action, under the explicit class and fibre hypotheses below.
For overlapping polynomial-direction outputs this puts the generic
image equation in the multiplier variables alone.

Claims HC4MYIP0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). The results concern
arbitrary degrees and all admitted stable dimensions. Both plane
endpoints and unrestricted HC4/JC2 remain open.

## 1. The two-sided component argument

Let `S` be a domain over an algebraically closed field `Omega`
of characteristic zero, and let `D` be a locally nilpotent
`Omega`-derivation. Define, for nonzero `a`,

```
delta(a)=max{n : D^n(a)!=0}.
```

Then `delta(ab)=delta(a)+delta(b)`, `delta(a+b)<=max(delta(a),delta(b))`
when the sum is nonzero, and `delta(c)=0` for `c in Omega*`.
For the product assertion, if the two degrees are `m,n`, the
Leibniz formula gives

```
D^(m+n)(ab)=binomial(m+n,m) D^m(a) D^n(b)!=0,
```

and the next derivative vanishes. This proves additivity directly.
In particular all units have degree zero. Degree zero is precisely
`ker D`, which is not assumed to be a field.

These are standard facts; [Alhajjar, arXiv:1310.3784v1,
Sections 1--2](https://arxiv.org/html/1310.3784v1) supplies primary
context and the usual derivation-filtration terminology. Their
elementary proof is reproduced here. No classification of locally
nilpotent derivations, finite-generation theorem for their kernels,
or Makar-Limanov computation is imported.

Suppose `e_i,e_j in S` generate distinct disjoint principal prime
divisors, with

```
(S/(e_i))*=(S/(e_j))*=Omega*,
delta(e_i)=delta(e_j)=m>0.                         (1)
```

Disjointness and the quotient unit groups give

```
e_j-c_ij=b_ij e_i, e_i mod e_j=c_ji,
c_ij,c_ji in Omega*.
```

Degree additivity gives `delta(b_ij)=0`. Reducing the first
identity modulo `e_j` gives

```
b_ij mod e_j = -c_ij/c_ji = b0 in Omega*.
```

Thus `b_ij-b0` belongs to `(e_j)`. If it were nonzero, its
degree would be zero, while any nonzero multiple of `e_j` has
degree at least `m>0`. This is impossible. Therefore `b_ij=b0`
is scalar, and `e_i,e_j` are affine-linear combinations of one
another. After a common coefficient normalization their nonconstant
parts agree.

This uses both affine-space components, not just the first
disjointness identity. It is exactly what allows a large
degree-zero kernel instead of the old scalar degree-zero ring.

## 2. One moving additive action gives the carried inverse (HC4MYIP0)

Let `K:A2->A2` be complex plane Keller and let a genuine closed
polynomial embedding `j:A^(r+2)->A^(r+3)`, `r>=0`, carry it
through ambient observables `Pi`. Put `B=C[T,V]`,
`F=Frac(B)`, `R=C[ambient variables]`, and let `E` be the
prime image equation. Assume:

1. `R` is faithfully flat over `B` and every affine height-one
   generic target fibre is geometrically integral.
2. `R_Omega`, `Omega=F_bar`, is a normal domain and its geometric generic image components are
   principal. Zero source-fixed class subgroup suffices for
   principality, as in HC4MYIN0.
3. There is an `F`-locally nilpotent derivation `D` of `R_F`
   with `D(E)!=0`.

Then `K` is a polynomial automorphism. No generic positive
filtration, global free filtration, slice, free additive action,
or finite-dimensional degree pieces are required.

Over `Omega`, the image equation factors as a unit times
distinct principal component generators `e_1,...,e_d`, with
`d=[C(s,t):F]`. Each quotient is a polynomial algebra in the
stable variables, so its unit group is `Omega*`. The components
are conjugate. A component generator may initially change by a
nonconstant unit, but every such unit has derivation degree zero.
Thus its degree is independent of that ambiguity. Since `D` is
defined over `F`, all conjugate component degrees agree. The
residual unit also has degree zero, so `D(E)!=0` makes their
common degree positive.

In fact this **forces** all geometric generic units to be scalar.
For any unit `w`, restriction to one component is `c in Omega*`.
If `w-c` were nonzero, its degree would be zero, while its
membership in `(e_i)` would give strictly positive degree.
Therefore `w=c`. No hypothesis on the generic unit group is needed.

Now choose
an arbitrary `F`-basis of `R_F` containing `1` and normalize
each generator's common last nonzero coefficient to one. This
is possible without a filtration: each support is finite, conjugacy
preserves its last position, and ambiguity of a generator is
only by a scalar.

The normalized generators still have equal positive degree. Section 1 now
gives `e_i=h0-alpha_i`, with `h0 in R_F` by Galois descent.
Its positive derivation degree makes it transcendental over `F`.

The rank-two reflexive pencil descent in
[HC4MYIN0](REFLEXIVE_PENCIL_DESCENT.md#1-a-primitive-global-coordinate-for-any-generic-pencil)
needs no filtration once `h0` has been supplied. It produces
`h in R` with `F[h]=F[h0]` and `F[h] intersection R=B[h]`.
Hence the literal global image equation is

```
E=B0(T,V,h), B0 in B[S], deg_S B0=d.                (2)
```

The substitution is injective. A factor of positive `S` degree
has positive derivation degree after substitution and cannot be
a unit; a nonconstant factor from `B` stays nonconstant in the
polynomial ambient ring. Thus `B0` is irreducible. On the image,
`h` is a constant label on every geometric generic stable fibre,
so `h composed with j=a(s,t)` is independent of the stable variables.
It has minimal relation `B0` over `F`.

The retained smooth-image and Keller-differential argument makes
`(B0)_S(K,a)` nowhere zero, hence a nonzero complex constant `c`.
The lower-label-degree relation `(B0)_S-c` then vanishes identically,
forcing `deg_S B0=1`. The birational Keller endpoint gives the
polynomial inverse. This proves the assertion in all degrees.

## 3. Trivial relative invariants close the output; otherwise the image is invariant (HC4MYIP1)

Under conditions 1--2 of Section 2, if `K` is nonautomorphic then

```
D(E)=0 for every D in LND_F(R_F).                   (3)
```

This is the direct contrapositive of HC4MYIP0, not an assumption
that arbitrary source changes preserve a Hessian. In particular
the image is invariant under every such relative additive action.
For an action defined only on `R_F`, multiply its derivation by
one nonzero element of `B` clearing the coefficients on finitely
many ambient generators. The multiplier is in its kernel, so
the resulting derivation of `R` is still locally nilpotent. Thus
the generic invariance also gives actual global relative additive
actions after this rescaling. It may degenerate at parameter levels;
no freeness or linearization there is inferred.

Assume condition 1, a geometric generic normal domain and principal geometric
generic image components. If
an explicitly supplied family of `F`-locally nilpotent derivations
has joint kernel exactly `F` on `R_F`, then `K` is invertible.
A finite supplied family suffices; no classification of all actions
is required. In particular a trivial relative Makar-Limanov
invariant is sufficient, where this invariant means
`intersection_(D in LND_F(R_F)) ker D`.

To check scalar geometric units, extend the supplied derivations
to `Omega`. Joint kernels commute here with scalar extension:
write any element using finitely many `F`-linearly independent
coefficients of `Omega`. Its annihilation by every derivation
forces every `R_F` coefficient into the original joint kernel.
Hence the joint kernel after extension is `Omega`. Every unit
is killed by every locally nilpotent derivation by degree additivity,
so all geometric generic units are scalar.

The image equation is not in `F`: its generic zero set is nonempty.
At least one supplied derivation therefore moves it. HC4MYIP0
applies. This replaces the generic-filtration requirement by a
joint-kernel certificate; neither condition is asserted equivalent
to the other, or automatic for arbitrary observables.

## 4. Every surviving overlapping polynomial-direction image is vertical (HC4MYIP2)

Consider

```
Pi=(T,M(T,X)Y-P(T,X,Z)),
p=#X>=1, q=#Z>=1, r>=p+q-1,
```

with unused variables `W`. Assume the HC4MYIO1 residual-direction
condition: modulo every geometric factor of `M` over `C(T)_bar`,
`P` retains positive joint `Z` degree. Assume also the exact
HC4MYIO2 content-level condition: at every zero `c` of the gcd
of the `X` coefficients of `M`, the full polynomial `P(c,X,Z)`
is nonconstant closed.

These hypotheses give generic smoothness, zero source-fixed class
subgroup, faithful flatness and geometrically integral height-one
fibres, without assuming any generic filtration. They also give
scalar geometric generic units directly, although HC4MYIP0 no longer
requires that as an input. In the localization
`Omega[X,1/M,Z,W]`, any unit of the unlocalized ring is
`c product a_j^(n_j)`. The exact row-divisor formulas
`div(a_j)=sum_i D_(j,i)`, with nonempty disjoint rows, force
every exponent to be zero. Hence the unit is scalar.

For each `Z_i` there is a global locally nilpotent derivation

```
D_i(T)=D_i(X)=D_i(W)=0,
D_i(Z_j)=delta_ij M,
D_i(Y)=partial P/partial Z_i.                       (4)
```

It preserves both observables. Its local nilpotence follows from
`D_i^n(Y)=M^(n-1) partial^n P/partial Z_i^n` for `n>=1`.
Include the ordinary derivatives in the unused `W` variables.
After inverting `M`, these become `M partial/partial Z_i`
and `partial/partial W_j`, with joint kernel `F[X,1/M]`.
Their joint kernel in the original generic ring is exactly

```
R_F intersection F[X,1/M] = F[X].                  (5)
```

Indeed every removed divisor dominates its multiplier component
by HC4MYIO1. A numerator not divisible by a multiplier factor
has valuation zero on its row divisors, while that factor has
valuation one. A denominator using that factor cannot be regular.
This removes every possible denominator in `M`, proving (5).
The same argument after geometric extension, followed by descent,
handles factors not defined over `F`.

If a carried `K` is nonautomorphic, HC4MYIP1 forces every
derivation (4) and every unused-variable derivative to kill `E`.
Therefore

```
E in F[X].                                        (6)
```

This is a uniform generic normal form for **all** survivors of
the overlap theorem, including those with no verified positive
generic filtration. Their image equation depends on `Y,Z,W`
only through the target observable `V`. A candidate whose actual
image equation violates (6) already forces the Keller inverse.

### 4.1 A global equation for primitive multipliers

If the coefficient content of `M` is constant, then

```
F[X] intersection R = C[T,V,X].                    (7)
```

At every height-one target prime `pi`, `M` stays a nonzero
polynomial in `X`: this is immediate horizontally and follows
from primitive content at `T=c`. Localizing the fibre at `M`
eliminates `Y`, so the map from `kappa(pi)[X]` to that fibre
is injective. For a polynomial in `F[X]` with a negative minimum
coefficient valuation, compensate its denominator at `pi`; its
nonzero reduction in `kappa(pi)[X]` cannot vanish in the fibre.
It therefore cannot have been integral in `R_pi`. Repeat at
every prime to prove (7).

Thus in the primitive case the actual polynomial image equation
has the global form

```
E(T,M Y-P,X) with E(T,V,X) in C[T,V,X].              (8)
```

No global free basis, polynomial quotient of the stable variables,
or rectification of the image is asserted by (8). For nonprimitive
content the generic conclusion (6) remains valid; descent (7)
requires another injectivity check at its content levels.

### 4.2 What the induced source actions fix

Under nonautomorphy, (4) annihilates the actual global image equation,
so descends to the polynomial source ring `C[s,t,w]`. It fixes
`F=C(f,g)`. Since `C(s,t)/F` is finite separable, any derivation
fixing `F` kills `s,t` themselves: differentiate their separable
minimal equations. Thus all these induced additive actions move
only stable variables. This is an actual-source constraint, with
no assumption that the actions are free, linear or simultaneously
rectifiable.

## 5. Controls and remaining frontier

The degree-zero hypothesis really has weakened. In `F[x,y]`,
`D=partial/partial y` has degree-zero ring `F[x]`; its degree-zero
piece is infinite-dimensional over `F`. Nevertheless conjugate
components `y-alpha_i` satisfy Section 1 and form the required
pencil. A degree-zero image equation such as `x` is invariant
under this action; one action alone does not exclude it.
Supplying also `partial/partial x` makes the joint kernel `F`.

The two-sided condition also matters. In `Omega[x,y]`, the
disjoint prime divisors `x=0` and `xy-1=0` have
`xy-1=-1+y*x`, with nonconstant multiplier `y`. The second
quotient is a Laurent polynomial ring and has nonconstant units,
so it fails Section 1. No affine-space-component conclusion is
being extended to this control.

For a suspension control, `M=X^2-T` and
`P=(X+1)Z^3+Z^2+X Z` have the derivation (4) in every degree.
Its generic kernel is `F[X]`; retaining that helper variable
is exactly the remaining invariant-image possibility. No statement
that the complete relative Makar-Limanov invariant equals `F[X]`
is made, since additional actions could reduce it further.

The global-intersection qualification has a sharp control:
`M=T`, `P=X`. Its content level `P(0,X,Z)=X` is closed and
the generic fibre is polynomial. Yet
`(X+V)/T=Y` belongs to `R intersection F[X]` and is not in
`C[T,V,X]`. Thus the primitive-content condition in (7) cannot
simply be dropped. This example satisfies the content condition;
it is not a Keller counterexample.

For an actual HC4 compatibility check retain
`M=1+X1 X2`, `P=X1 Z^m+Z+X1 C(T,X2,Z)` with `deg_Z C<m`
and the closed identity embedding from HC4MYIO2. Formula (4)
preserves both observables, and its derivation degree on the
identity-image equation `X1` is zero. This illustrates why a
stationary image under one action does not imply nonautomorphy.
The actual unit-Hessian quadratic pullback remains unchanged.

Remaining work concerns additive-invariant images (6), whether
other relative actions or actual-source geometry close them,
nonconstant composite content levels, other observable geometries,
compatible observables for arbitrary HC4 graphs, actual endpoint
representatives, both plane endpoints and unrestricted HC4/JC2.
The former lack-of-generic-filtration residual has been reduced
to a specific invariant image geometry under the stated hypotheses;
it has not been declared universally absent.

The [source record](sources/additive-action-context.json),
[checker](verify_additive_actions.py),
[receipt](receipt-additive-actions.json) and
[validation](validation-additive-actions.json) preserve the exact
scope and prior proof bytes. Formal verification, independent replay,
external review and novelty priority are not claimed.
