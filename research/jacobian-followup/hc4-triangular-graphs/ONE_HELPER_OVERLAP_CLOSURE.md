# Every one-helper polynomial-direction overlap closes

15 November 2026. The additive-invariant image reduction closes
completely when there is one multiplier variable. The rank-two
descent lemma applies to that variable itself; there is no need
for its original coefficients to descend globally, or for the
map forgetting root variables to be flat. The smooth-image
derivative forces the source helper into the Keller target field.
Generic elimination, including the multiplier-zero case, then
forces degree one for the original plane field extension.

Claims HC4MYIQ0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). All polynomial,
coefficient, embedding and stable degrees are unrestricted. Both
plane endpoints and unrestricted HC4/JC2 remain open.

## 1. A generic univariate image equation has a rational source helper (HC4MYIQ0)

Let a genuine polynomial closed embedding
`j:A^(r+2)->A^(r+3)`, `r>=0`, carry a complex plane Keller
pair `K=(f,g)` through ambient observables `Pi=(T,V)`. Put
`B=C[T,V]`, `F=Frac(B)`, `R=C[ambient variables]`, and
`L=C(s,t)`. Let `E` be the prime image equation, so the
actual source ring is `R/(E)=C[s,t,w_1,...,w_r]`.

Assume `R` is faithfully flat over `B` and every affine
height-one generic target fibre is geometrically integral. Suppose
there is an element `h0 in R_F` transcendental over `F` with

```
E=Q(h0), Q in F[S].                                (1)
```

Then `Q` has degree one, and the restriction of `h0` to the
generic source is in `F`. No generic positive filtration, principal
geometric component hypothesis, or additive action is needed once
(1) is supplied. This statement concerns the helper label; its
degree one does not by itself assert `[L:F]=1`.

### 1.1 Descend the helper before differentiating

The rank-two lattice lemma in
[HC4MYIN0](REFLEXIVE_PENCIL_DESCENT.md#1-a-primitive-global-coordinate-for-any-generic-pencil)
applies to `h0`, independently of how the pencil was obtained.
It provides `h in R` with

```
F+F h=F+F h0,
F[h] intersection R=B[h].                          (2)
```

Thus the actual global image equation has the literal expression

```
E=C(T,V,h), C in B[S], deg_S C=deg Q>0.             (3)
```

The degree is positive because a nonzero element of `F` is a
generic unit, whereas the geometric generic image is nonempty.
The substitution `B[S]->R` is injective since `h` is
transcendental over `F`. Every nonconstant factor of `C` therefore
substitutes to a nonconstant element of the polynomial ring `R`,
whose only units are `C*`. Irreducibility of `E` makes `C`
irreducible in `B[S]`, hence in `F[S]` by Gauss.

On the source, `h` satisfies (3), so it is algebraic over `F`.
It belongs to the polynomial ring `L[w]`, whose nonconstant
polynomials in `w` are transcendental over `L`, and therefore
cannot be algebraic over `F`. Hence

```
h composed with j=a(s,t) in C[s,t].                 (4)
```

The irreducible polynomial `C` is its minimal relation over `F`.

### 1.2 The smooth derivative makes the label linear

If `C_S(K,a)` vanished at a source point, differentiating (3)
and pulling back to the source would give
`C_T df+C_V dg=0` there. The invertible Keller differential
forces `C_T=C_V=0`. All three coefficients in the ambient
identity `dE=C_T dT+C_V dV+C_S dh` would then vanish at the
image point, contradicting smoothness of the actual hypersurface
image. Thus `C_S(K,a)` has no zeros and is a nonzero complex
constant `c`.

The lower-label-degree polynomial `C_S-c` vanishes at `a`.
Minimality gives `C_S=c` identically. Characteristic zero yields
`C=c S+b(T,V)`. Thus `deg Q=1`, `a=-b(K)/c`, and the
affine `F`-relation between `h0` and `h` puts `h0` itself in
`F` on the generic source. This proves the claim. It is the
retained smooth minimal-label mechanism, now applied after
descent of a supplied univariate image equation rather than
construction of a whole generic component pencil.

In particular, under these hypotheses an actual image equation
cannot be a nontrivial polynomial composition `Q(h0)` with
`deg Q>=2`. Special geometric components of its zero fibre
may still require a separate field argument; that is supplied
below for one-helper overlap outputs.

## 2. All one-multiplier-variable polynomial-direction outputs (HC4MYIQ1)

Take one multiplier variable `X`, any `q>=1` root variables
`Z=(Z_1,...,Z_q)`, and arbitrary polynomials

```
M(T,X) in C[T,X], M!=0,
P(T,X,Z) in C[T,X,Z].
```

Require precisely the retained polynomial-direction and content
conditions:

1. Over `k_bar`, `k=C(T)`, reduction of `P` modulo every
   geometric irreducible factor of `M` has positive joint `Z`
   degree. If `M` is a generic parameter unit, there are no rows.
2. At every root `c` of the gcd `delta_M(T)` of all `X`
   coefficients of `M`, the full polynomial `P(c,X,Z)` is
   nonconstant and closed over `C`.

Then every genuine closed stable embedding carrying a plane
Keller pair through

```
Pi=(T,M(T,X)Y-P(T,X,Z)), r>=q,                     (5)
```

forces its polynomial inverse. Include `r-q` unused variables.
All coefficient, parameter, multiplier, root, source and embedding
degrees, root-variable counts and admitted stable dimensions are
unrestricted. There is no generic-filtration, monicity, Bezout
leading, primitive-content, special source-fibre or forgetful-flatness
hypothesis. In particular the admitted nonprimitive content level
may depend only on `X`, with no surviving `Z` term there.

### 2.1 The helper becomes rational over the target

Suppose for contradiction that `K` is nonautomorphic. Conditions
1--2 give the full hypotheses of HC4MYIP2. It puts the actual
image equation, over the generic target field, in `F[X]`.
The ambient `X` is transcendental over `F`: after inverting the
nonzero multiplier, the generic ring is `F[X,1/M,Z,W]`.

Apply HC4MYIQ0 with `h0=X`. It proves that the image equation
has degree one in `X` over `F` and that the actual source
restriction is

```
X=a in F=C(T,V).                                   (6)
```

This uses the rank-two saturation even when the old global
intersection `F[X] intersection R=B[X]` fails at a content level.
We do not need to assert that the original equation was already
polynomial over `B[X]`.

As an ideal of the generic ring the image equation is therefore
`(X-a)`. Its quotient is also, from the actual embedding,

```
R_F/(E)=L[w_1,...,w_r].                            (7)
```

Indeed localizing `C[s,t]` by the nonzero target polynomials
gives a finite-type algebraic `F`-domain, hence the field `L`.

### 2.2 The multiplier does not vanish on the image

If `M(T,a)!=0` in `F`, it is a field unit. Substitute (6)
in the generic ambient equation and eliminate `Y`, giving

```
R_F/(E)=F[Z,W].                                   (8)
```

The field `F` is relatively algebraically closed in a purely
transcendental extension. Comparing (7)--(8) puts the algebraic
field `L` inside `F`, hence `L=F`. The retained birational
Keller case gives the polynomial inverse, a contradiction.

### 2.3 The multiplier vanishes on the image

If `M(T,a)=0`, then `a in k(V)` is algebraic over `k`, so
`a in k=C(T)`. To see the latter directly, a nonconstant
rational function of an indeterminate `V` is transcendental
over `k`; its rational equation would otherwise make `V`
algebraic over `k`. Define

```
P_a(T,Z)=P(T,a(T),Z) in C(T)[Z].
```

Condition 1 ensures it is nonconstant. The generic image is now

```
R_F/(E)=F[Y,Z,W]/(P_a+V).                          (9)
```

Base change (9) from `F` to the original source field `L`.
The multivariable polynomial-output irreducibility theorem
HC4MYIL0, also included in HC4MYIO0, makes

```
L[Y,Z,W]/(P_a(f,Z)+g)
```

a domain. But base changing (7) gives
`(L tensor_F L)[w]`. A nontrivial finite extension cannot have
`L tensor_F L` a domain: as a finite-dimensional `L`-algebra,
a domain would be a field, and its surjective multiplication
map to `L` would be injective. Its dimension over `L`, equal
to `[L:F]`, would therefore be one. This proves `L=F` here
as well. No Galois assumption on `L/F` is used.

For primary context, the standard separable self-tensor description
is in [Stacks 05DS, Lemma 10.43.8 and its proof](https://stacks.math.columbia.edu/tag/05DS).
The simpler finite-dimensional-domain argument needed here is
fully reproduced above and does not import that lemma.

Both cases contradict nonautomorphy. This proves (5), including
the multiplier-zero case, arbitrary nonprimitive content admitted
by condition 2, and every root-variable count. Ambient and target
polynomial automorphism orbits inherit the carried-output result
by the retained change-of-observables argument. No nonlinear
source Hessian invariance is used.

## 3. The remaining HC4 overlap has two helpers and one root variable (HC4MYIQ2)

Within the polynomial-direction overlap family of HC4MYIO1/2,
HC4 size means `r=2`. Its variable count is

```
p=#X>=1, q=#Z>=1, u=#W>=0,
p+q+u=3.
```

Every `p=1` case closes by HC4MYIQ1. Thus any nonautomorphic
carried HC4 pair in this admitted family must have

```
p=2, q=1, u=0.                                    (10)
```

HC4MYIP2 still forces its image equation into `F[X_1,X_2]`,
with a literal global equation in `C[T,V,X_1,X_2]` when the
multiplier content is primitive. HC4MYIQ0 further excludes any
nontrivial polynomial composition of the image equation in a
transcendental generic ambient helper. This is a remaining
invariant plane-curve geometry, not a resolved inverse theorem.

Condition (10) is a structural variable count and imposes no
degree bound. It is not a universal normal form for arbitrary
HC4 graph embeddings. Residual-direction failures, nonconstant
composite content levels, other output geometries, compatible
arbitrary-graph observables and actual endpoint representatives
remain separate. Both plane endpoints and unrestricted HC4/JC2
remain open.

## 4. Exact controls and boundaries

A uniform admitted overlap, with different degrees on different
multiplier components, is

```
M=X(X-1), P=(1-X) Z^m+X Z^n+T X Z,
m,n>=1.
```

At `X=0` and `X=1` the residual degrees are respectively
`m,n`; no constant leading coefficient common to both rows is
required. The theorem covers all `m,n` at once and does not
assume a positive generic filtration has been constructed.

A nonprimitive content control is

```
M=T X, P=X+T Z^m, m>=1.
```

The generic multiplier row `X=0` retains `T Z^m`, and at its
content zero `T=0`, `P=X` is nonconstant closed although
independent of `Z`. Here `(X+V)/T=Y X-Z^m` is a genuine
element of `R intersection F[X]` outside `B[X]`. The rank-two
coordinate in Section 1 is precisely what avoids that old
coefficient-denominator obstruction.

The two generic cases in Sections 2.2--2.3 must not be conflated.
When `M(T,a)=0`, `Y` does not eliminate. That case needs the
retained polynomial-output irreducibility and the self-tensor
argument; a bare assertion that a degree-one helper makes the
entire original source extension degree one would be invalid.

Likewise a rational source helper alone is not global polynomial
reconstruction, and does not invoke Formanek's criterion by
itself. The proof uses the explicit generic ambient equation to
finish the field argument. It does not assume a helper quotient
is a polynomial ring or use stable cancellation.

A non-Keller control shows the role of the smooth derivative:
the closed embedding `(s,t,w)->(s,t,0,w)` through
`Pi=(T,X^2-Y)` carries `(s,t^2)`. Its image equation `Y`
becomes `X^2-V` over `F[X]`, of helper degree two. Its source
Jacobian vanishes at `t=0`, exactly where the label derivative
`2t` vanishes. No Keller conclusion applies.

A genuine HC4 control uses `M=1+(X-1)^2` and
`P=T Z_1^3+X Z_2^2+Z_1 Z_2`, with closed embedding

```
j(s,t,u,v)=(T,X,Y,Z_1,Z_2)
         =(s,1,t+P(s,1,u,v),u,v).
```

Its inverse on the image is `(T,Y-P(T,X,Z),Z_1,Z_2)`.
It carries the identity pair. The independent nondegenerate
ambient quadratic `T Z_2+Y Z_1+X^2` pulls back to
`1+s v+t u+u P(s,1,u,v)` and has actual Hessian determinant
one. This is compatibility with a closed embedding and a
unit-Hessian potential, not a nonautomorphic realization.

The [source context](sources/one-helper-overlap-context.json),
[checker](verify_one_helper_overlap.py),
[receipt](receipt-one-helper-overlap.json) and
[validation](validation-one-helper-overlap.json) preserve the
written general proof separately from exact controls and byte
pins. No new unproved external theorem, action classification,
forgetful-flatness assertion, cancellation theorem, formal
verification, independent replay, external review or novelty
priority is claimed.
