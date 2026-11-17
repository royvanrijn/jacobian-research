# A failed residual direction can retain the entire plane Keller problem

15 November 2026. The residual polynomial-direction condition cannot
simply be discarded from the carried-inverse theorem. Every plane
Keller pair has a canonical closed HC4-size embedding through an
observable map with smooth geometrically integral fibres over the
whole target. Its geometric generic class group records the original
field degree exactly. The same embedding has actual unit Hessian for
a nondegenerate ambient quadratic, and transfers to the fixed
Meng--Yang potential through the retained ambient normal form.

This is a sharpness and equivalence result, not a proof of JC2 or a
nonautomorphic example. The source potential constructed below has
an explicit polynomial gradient inverse for every input pair. The
carried pair must not be identified with the intrinsic plane pair of
a hypothetical nonautomorphic gradient.

Claims HC4MYIS0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). Both plane endpoints,
compatible observables for arbitrary intrinsic HC4 data, and
unrestricted HC4/JC2 remain open.

## 1. A universal closed graph with good target fibres (HC4MYIS0)

Let `K=(f,g):A2->A2` be any complex plane Keller pair, with
`J=f_1 g_2-f_2 g_1 in C*`. On the polynomial ambient five-space
with coordinates `(T,X1,X2,Y,Z)`, set

```
M=f(X1,X2)-T,
V=M Y+g(X1,X2),
Pi=(T,V),
j(s,t,u,v)=(f(s,t),s,t,u,v).                        (1)
```

This is an actual closed embedding, with inverse `(X1,X2,Y,Z)`
on its image and smooth prime image equation `E=M`. It carries
`K(s,t)` through `Pi`. In the overlap notation it has two
multiplier helpers, one unused root variable `Z`, and
`P=-g(X1,X2)`. The multiplier coefficient content is primitive.
The residual polynomial direction fails on the entire multiplier
row, because `P` has no `Z` term at all.

Nevertheless the target fibres are better than the earlier sufficient
height-one condition. For any algebraically closed extension `k/C`
and any `a,b in k`, the fibre is

```
k[X1,X2,Y,Z]/((f(X)-a)Y+g(X)-b).                   (2)
```

The two coefficients `f-a,g-b` are coprime in `k[X1,X2]`.
If an irreducible `h` divided both, reduction modulo `h` would make
their two derivative rows multiples of `dh`, so their determinant
would vanish modulo `h`, contradicting `J in C*`. Thus the displayed
linear polynomial in `Y` is primitive and irreducible by Gauss's
lemma. Its quotient is a domain of dimension three. It is nonempty:
choose any `X` with `f(X)!=a` and solve for `Y`, with arbitrary `Z`.

Put `H=(f-a)Y+g-b`. Its two helper derivatives satisfy the
literal polynomial identity

```
f_1 H_2-f_2 H_1=J.                                (3)
```

Hence their ideal is the unit ideal and every geometric fibre (2)
is smooth. The argument also works over the algebraic closure of
the residue field of any target point, not just closed complex
points. All target fibres are geometrically integral and smooth,
with their expected dimension. In particular the morphism is
faithfully flat: the regular-base, Cohen--Macaulay-source dimension
criterion is the retained
[Stacks 00R4 miracle-flatness lemma](https://stacks.math.columbia.edu/tag/00R4).
Surjectivity was checked directly. The geometric generic ambient
ring is a smooth normal domain.

Thus faithful flatness, even smooth geometrically integral fibres
everywhere, does not by itself replace the missing residual
polynomial direction in the class calculation.

## 2. Its generic class group retains precisely the source degree (HC4MYIS1)

Put `B=C[T,V]`, `F=C(T,V)`, `L=C(s,t)` via `T=f,V=g`,
`d=[L:F]`, and let `Omega` be an algebraic closure of `F`
containing this copy of `L`. The geometric generic ambient ring is

```
S=Omega[X1,X2,Y,Z]/(M Y+g-V).                      (4)
```

The retained regular-parameter theorem HC4MYHX0 says that
`C(f)` is relatively algebraically closed in `C(s,t)`. Consequently
`f(X)-T` is absolutely irreducible over `C(T)` and remains
irreducible over `Omega`. Here this is the generic fibre over `T`;
it is not a claim that all special fibres of `f` are irreducible.

Inverting `M` gives exactly

```
S[1/M]=Omega[X1,X2,1/M,Z],                         (5)
```

a UFD whose units are `Omega* M^Z`. Without inverting `M`,

```
S/(M)=Omega[X1,X2]/(f-T,g-V)[Y,Z]
      =product_(i=1)^d Omega[Y,Z].                 (6)
```

Indeed before extending `F`, the helper quotient is exactly the
finite separable field `L`: the generic Keller source algebra is an
algebraic finitely generated `F`-domain and hence a field. There
are precisely `d` reduced geometric points. Thus `M=0` in `S`
has `d` disjoint prime divisors `D_i`, each the actual source
component `A2`. Reducedness of (6) gives valuation one at every
such divisor, so

```
div_S(M)=D_1+...+D_d.                              (7)
```

The elementary localization argument now gives the complete lattice

```
Cl(S)=Z^d/Z(1,...,1).                              (8)
```

For completeness, every Weil class has a representative supported
on the `D_i`, since (5) is factorial: subtract a rational principal
divisor whose restriction represents its divisor on that open set.
Any principal relation supported there becomes a unit in (5).
Those units are exactly `c M^n`, giving only the multiples of (7).
This proves both generation and all relations, rather than only
a rank bound. The ring is normal, so a rational function with zero
divisor on an affine open is a regular unit there. The same reasoning
shows

```
S*=Omega*:                                        (9)
```

a global unit is `c M^n` after localization, and (7) forces `n=0`.
No general class-group classification is imported.

Let `H=Gal(Omega/L)` and select the source component `D_1`
corresponding to the fixed embedding `L->Omega`. It is fixed by
`H`, so its class is source-fixed. If `d>1`, its vector `e_1`
is nonzero of infinite order modulo the diagonal in (8). Therefore

```
Cl(S)^H=0  <=>  d=1  <=>  K is a polynomial automorphism.
                                                               (10)
```

The last equivalence uses the retained birational Keller theorem.
More precisely, if `m` is the number of `H`-orbits on the `d`
embeddings, then `Cl(S)^H` is free of rank `m-1`. If a class
represented by an integer vector is fixed, its difference from each
permuted vector is an integer multiple of the diagonal. Summing
coordinates makes that multiple zero. Thus its representative is
constant on each orbit; quotienting by the diagonal gives the rank
claim. When `d>1` the selected singleton orbit leaves at least one
other orbit. No Galois hypothesis on `L/F` is used.

The global equation `M` is principal and smooth, but its individual
geometric generic components are nonprincipal when `d>1`.
The global equation being linear in one helper label is not enough
to make the original source field degree one.

The map forgetting `Y,Z` also exhibits the precise failure of the
earlier two-helper proof. Off `M=0` its fibre is `A1_Z`. On the
helper curve `M=0` it has empty fibres except at the `d` points
`f=T,g=V`, where it has full fibres `A2_(Y,Z)`. It is not flat.
The removed divisors project to these points, not dominantly onto
the multiplier curve. The scalar-unit result (9) survives, but the
IO1 row-dominance and source-transitivity conclusions do not.

## 3. Unit Hessian does not remove this sharpness boundary (HC4MYIS2)

Take the nondegenerate ambient quadratic

```
Q=T^2/2+X1 Y+X2 Z.
```

Its pullback under the literal parameterization (1) is

```
phi(s,t,u,v)=f(s,t)^2/2+s u+t v.                   (11)
```

For any polynomial `A(s,t)`, the Hessian of `A+s u+t v` has
block form `[[Hess A,I],[I,0]]` and determinant one. In this case
the gradient and its polynomial inverse are explicit:

```
grad phi=(u+f f_s, v+f f_t, s, t),
(a,b,c,d0) -> (c,d0,a-f(c,d0)f_s(c,d0),
                         b-f(c,d0)f_t(c,d0)).     (12)
```

This holds in every degree and does not use invertibility of `K`.
In particular it supplies no nonautomorphic HC4 gradient.

These are also genuine closed embeddings for the fixed scaled
Meng--Yang potential. HC4MYHM0 supplies the retained polynomial
ambient automorphism `Phi` with
`Psi=lambda_*+Q5 composed with Phi`, where
`Q5=L T0^2+N y0 Z0+2N x0 W0`, `LN!=0`.
Choose a constant `c` with `L c^2=1/2` and the invertible linear
map

```
ell(T,X1,X2,Y,Z)=(x0=X2,y0=X1,T0=cT,
                                      Z0=Y/N,W0=Z/(2N)).
```

Then `Q5 composed with ell=Q`. The embedding
`Phi^(-1) composed with ell composed with j` is closed, has
the same actual source potential `lambda_*+phi`, and carries
`K` through `Pi composed with ell^(-1) composed with Phi`.
This is a literal function identity in the original source
parameters. No nonlinear source Hessian invariance is assumed.

Consequently each of the following universal assertions is
equivalent to plane JC2:

1. Every carried pair from the canonical family (1) is invertible.
2. Every two-helper overlap with primitive multiplier content and
   smooth geometrically integral target fibres forces the carried
   inverse even when the residual polynomial direction is dropped.
3. The assertion in item 2 restricted to genuine HC4-size unit-Hessian
   embeddings, including ambient polynomial-automorphism orbits.
4. The canonical family (4) always has zero source-fixed class group.

JC2 immediately implies all four. Conversely (1) realizes every
plane Keller pair, with all the stated fibre and actual embedding
properties, and (10)--(12) give the other implications. These
are equivalences of open universal assertions, not assertions that
the inverse conclusions hold. The constructed potential always has
an automorphic gradient, so the carried pair need not be the
intrinsic coefficient pair of a nonautomorphic HC4 gradient.

## 4. Source context and the next useful gate

The generic surface underlying (4) is the affine modification

```
Omega[X1,X2][(V-g)/(f-T)],
```

followed by adjoining `Z`. This is standard terminology; primary
context is [Kaliman--Zaidenberg, arXiv:math/9801076v2,
Section 1](https://arxiv.org/html/math/9801076v2), specifically the
definition and the subalgebra `A[I/f]`. The displayed presentation,
all fibre claims, unit and class lattices, the source-fixed class
calculation and the unit-Hessian identities are proved above. No
topological preservation, exotic-space classification, flexibility
or cancellation theorem from that paper is imported.

The [source audit](sources/universal-modification-context.json),
[checker](verify_universal_modification.py),
[receipt](receipt-universal-modification.json) and
[validation](validation-universal-modification.json) separate
written all-degree reasoning from exact algebraic controls. No
nonautomorphic sample, formal verification, independent replay,
external review, or novelty priority is claimed.

A uniform extension across failed residual directions must overcome
this source-correlated class lattice, or exploit additional
compatibility with the intrinsic HC4 gradient data. Merely adding
smoothness, all-fibre integrality, scalar generic units, primitive
content, an actual closed embedding or actual unit Hessian cannot
remove the boundary. Compatible observables for arbitrary intrinsic
graphs, composite content levels, actual endpoint representatives,
both plane endpoints and unrestricted HC4/JC2 remain open.
