# Relative central kernels and all-degree binomial outputs

15 November 2026. The remaining quotient in the plane
[branch-inertia criterion](BRANCH_CURVES_AND_RATIONAL_TRIANGULAR_PARAMETERS.md)
need not be nilpotent. A relatively central kernel over any separately
preserved quotient suffices. For covers of a product with a punctured
parameter line, the quotient really is an algebraic parameter cover:
a section and full faithfulness of the complex finite-cover comparison
prove this, rather than an identification of abstract groups alone.

This closes every plane output `a(T)z^m+z^q+d(T)` in all degrees and,
more generally, `a(T)z^m+b(T)z^q+d(T)` whenever `b^m/a^q` has a finite
pole. The polynomial `a` may have arbitrarily many zeros and arbitrary
multiplicities; both `a` and `b` may vary. The complementary no-finite-pole
case is not asserted closed. Claims `HC4MYIJ0/1/2/3` are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json).

## 1. A relative central series over a preserved quotient (HC4MYIJ0)

Let `K:A^n->A^n` be complex Keller, `n>=2`, let `F=C(K)` and
`L=C(x_1,...,x_n)`, and let `M/F` be finite Galois with group `G`.
Write `H` for the image of `Gal(ML/L)` in `G`. Thus `M^H=M intersect L`.
Suppose a normal subgroup `N` is already certified to lie in `H`.
In `Q=G/N` supply a chain of normal subgroups

```
1=C_0 <= C_1 <= ... <= C_k=C,
C_i/C_(i-1) subset Z(Q/C_(i-1)).              (1)
```

Require that the source image preserve the full quotient `Q/C`.
Then `H=G`: the full group of `M` is preserved. All extension,
source and polynomial degrees are unrestricted.

Put `J=H/N`. Initially `JC_k=Q`. Whenever `JC_i=Q`, the subgroup
`JC_(i-1)/C_(i-1)` is normal in `Q/C_(i-1)`: every element of this
quotient is the product of an element of that subgroup and an element
of its centre, by (1). Thus `JC_(i-1)` is normal in `Q`. If it were
proper, its inverse image in `G` would be a proper normal subgroup
containing `H`. Its fixed field would be a nontrivial finite Galois
intermediate field inside `L`, forbidden by
[HC4MYIB0](NORMAL_QUOTIENT_AND_BOUNDARY.md#1-no-finite-galois-intermediate-field-hc4myib0).
It follows that `JC_(i-1)=Q`. Descending the entire finite chain gives
`J=Q` and `H=G`.

An ordinary central kernel is the length-one case. A nilpotent whole
group is another special case, with trivial final quotient. This
criterion permits an arbitrary nonnilpotent `Q/C` when that particular
quotient has its own preservation certificate. A nilpotent normal
kernel alone is insufficient: its conjugation by the quotient can
destroy (1). For example a transposition subgroup of `S3` surjects
onto `S3/A3` but is proper, and `A3` is a nilpotent normal kernel.
Here its nontrivial quotient action prevents relative centrality.

The argument is an all-dimensional field statement. A separate
geometric argument must supply `N subset H`; the plane-specific
nonparametrizable-curve test is one such source, not a higher-dimensional
consequence of (1).

## 2. Product covers with nilpotent fibre image (HC4MYIJ1)

Use polynomial target coordinates `T,X_2,...,X_n`. Let
`U=A1 minus S`, where `S` is any finite set, and let `Y` be a connected
Zariski open subset of `A^(n-1)`. Suppose the normalization of a finite
Galois `M/F`, `F=C(T,X_2,...,X_n)`, is finite etale over

```
W=U times Y.                                (2)
```

Let `C` be the image of the fibre factor `pi_1(Y)` in its finite
topological monodromy group `Q=Gal(M/F)`. If this **image** is nilpotent,
every complex Keller map in dimension `n` preserves the full group `Q`.
No nilpotence of the whole group is required. In particular this holds
for every such cover over

```
U times (C*)^h times A^(n-h-1), 0<=h<=n-1.    (3)
```

The open `W`, as well as its normalization cover, is connected and
nonsingular. The monodromy representation comes from
`pi_1(W)=pi_1(U) times pi_1(Y)`. The images of the two factors commute
and generate `Q`. Consequently `C` is normal, and the upper central
series of nilpotent `C` is a chain (1) in `Q`: the parameter-factor
image commutes with every element of `C`. In (3), `pi_1(Y)=Z^h`, so
`C` is even abelian and central.

### The quotient descends as an algebraic cover

Take the finite etale quotient cover `Z->W` with group `Q/C`, and
choose a complex point `y_0` of `Y`. Its restriction `Z_U` to the
algebraic section `T -> (T,y_0)` is an algebraic finite etale cover
of `U`. The quotient monodromy kills `pi_1(Y)`, so the restriction
is connected and its pullback `Z_U times Y` is isomorphic to `Z`
as a topological cover of `W`.

The required algebraic isomorphism is exactly **full faithfulness** in
[Milne, version 2.21, Theorem 3.4, printed page 28](https://www.jmilne.org/math/CourseNotes/LECc.pdf).
Both covers already exist algebraically. This is a new, explicitly
audited use of the theorem beyond the previously retained forward
covering assertion. Essential surjectivity, and the construction of
an algebraic cover from an arbitrary topological one, are not used.
The theorem applies because `W` is nonsingular over `C`.

It follows on function fields that

```
M^C = E(X_2,...,X_n),
E/C(T) finite Galois, Gal(E/C(T))=Q/C.        (4)
```

This is an equality inside the chosen algebraic closure, up to the
`F`-embedding supplied by the isomorphism; it is not only an abstract
group identification. Regularity of `C(f_1)` in the original Keller
source field, [HC4MYHX0](REGULAR_PARAMETER_MONODROMY.md), preserves
the entire finite parameter cover `E/C(T)`. Thus the source image
surjects onto `Q/C`. Apply Section 1 with `N=1`. This proves the claim.

All source and cover degrees, the size of `S` and the torus dimension
are unrestricted. Finite groups coming from `U` can be nonnilpotent;
they are not declared nilpotent because the fibre image is nilpotent.
The product structure and that fibre-image condition are actual
hypotheses. Removing arbitrary moving horizontal curves need not
leave a product, and its monodromy need not satisfy (1).

## 3. A plane branch quotient may have arbitrary parameter monodromy (HC4MYIJ2)

Now `K=(f,g)` is plane Keller, `F=C(T,V)` and `M/F` finite Galois.
Let `N_bad` be the normal subgroup generated by **actual divisorial
inertia** on target curves having no nonconstant polynomial `A1`
parametrization, as in
[HC4MYII3](BRANCH_CURVES_AND_RATIONAL_TRIANGULAR_PARAMETERS.md#4-nonparametrizable-branch-inertia-controls-all-factors-hc4myii3).
If the normalization in `M^N_bad` is finite etale over some

```
(A1 minus S) times (A1 minus {0}),            (5)
```

then the full group of `M/F` is preserved by `K`.

Indeed `N_bad subset H` by HC4MYII0/3. Section 2, in dimension two,
preserves the entire residual quotient cover, including its possibly
nonnilpotent parameter group. Combining these two facts gives `H=G`.

It is enough that every remaining finite branch curve lie in the
vertical lines `T in S` and one horizontal line `V=0`: the retained
purity theorem Stacks `0BMB` makes the finite normalization etale
over the regular complement (5). A polynomial target shear permits
the horizontal graph `V+d(T)=0` instead. Rational graph coordinates
are also allowed on a specified product open after removing their
poles from `U`, but the actual product certificate is still required.

The same reasoning applies to any residual quotient with the product
and nilpotent fibre-image certificate of Section 2. Full root-group
preservation gives irreducibility of any corresponding irreducible
separable polynomial. For filtered ambient outputs, equality of the
source and target class-action images permits the unchanged retained
[HC4MYHX1](REGULAR_PARAMETER_MONODROMY.md#2-nonabelian-monodromy-from-one-regular-target-parameter)
and HC4MYHZ2 global arguments. All filtration, normality, torsion-free
class-group and affine height-one fibre hypotheses remain required.

The assertion about `N_bad` is plane-only. Section 2 itself works
in every dimension. Neither section asserts full preservation of
arbitrary two-variable target covers.

## 4. A finite-pole test for two varying binomial coefficients (HC4MYIJ3)

Let

```
P(T,z)=a(T)z^m+b(T)z^q+d(T),
a,b nonzero polynomials, 1<=q<m, k=m-q.      (6)
```

If `b(T)^m/a(T)^q` has a finite pole, then every plane Keller pair
preserves irreducibility of `P(f,Z)+g` and its **full target splitting
group**. Every carried closed stable embedding through `(T,xy-P)`,
`r_stable>=1`, therefore forces the Keller inverse by HC4MYHZ2.
All degrees, leading multiplicities and numbers of leading levels
are unrestricted. The finite-pole condition is equivalently

```
there exists c with q ord_c(a)>m ord_c(b),
or equivalently a^q does not divide b^m.     (7)
```

After shearing the target coordinate, take `d=0`. At a nonzero
critical root, away from coefficient zeros,

```
m a z^k=-q b,
V=-k b z^q/m,
a(T)^q V^k = c_(m,q) b(T)^m,
c_(m,q)=(-1)^m q^q k^k/m^m !=0.             (8)
```

Reduce the last equation over `C[T]` by cancelling
`gcd(a^q,b^m)`. Every horizontal critical-value component is a
component of

```
U(T)V^k-W(T)=0, gcd(U,W)=1,
W/U=c_(m,q) b^m/a^q.                       (9)
```

Under (7), `U` is nonconstant. Bezout for `U,W` and (9) make
`U(T)` a unit in the coordinate ring of this entire curve, and
hence of every irreducible component. In a polynomial `A1` map,
`U(T(x))` is a polynomial unit, forcing `T(x)` constant. Then
`U(T(x))!=0` and (9) force `V(x)` constant as well. Thus none of
these components admits a nonconstant polynomial parametrization.
This argument works even when the critical-value equation is
reducible and the generic root action is imprimitive.

The root equation `P+V` is irreducible over `C(T,V)`: it is linear
and monic in `V` over `C[T,z]`, and is primitive over `C[T,V]`.
In characteristic zero it is separable. Outside the leading and
root-discriminant divisors the roots are integral and distinct, so
the finite splitting cover is etale there. The derivative is
`z^(q-1)(m a z^k+q b)`. Thus away from `a b=0`, every finite
branch lies on `V=0` or the critical-value curves (9). We do not
need to identify their inertia groups or classify the root group;
whichever actual branch components occur on (9) are all bad in the
sense of Section 3.

Kill `N_bad`. Every remaining branch curve is contained in
`a(T)b(T)=0` or `V=0`. Purity therefore gives (5), with `S` the
finite zero set of `a b`. Section 3 proves full splitting-group
preservation. Irreducibility and the global carried inverse follow.

### Consequences and exact boundary

In particular every

```
a(T)z^m+z^q+d(T), a!=0, 1<=q<m,             (10)
```

closes. Nonconstant `a` satisfies (7); constant `a` is already a
fixed-polynomial family, closed by HC4MYII2. This removes the
pure-power-leading restriction of the retained HC4MYII3 example.
It includes the previously closed adjacent case and all nonadjacent
exponents, with arbitrary leading coefficients.

Both coefficients may share zeros of arbitrarily large orders.
For example `a=T^(m+1)`, `b=T`, every `m>=2` and `1<=q<m`, has
`q(m+1)>m` and is admitted. No bounded degree, leading gcd, factor
degree, root primitivity or group classification is needed.

If `a^q` divides `b^m`, test (7) is unavailable: (9) can have a
polynomial parametrization. For example `a=T^m,b=T^q` gives
`P=(Tz)^m+(Tz)^q+d`, already closed by HC4MYII2 but not by the
finite-pole test. Thus failure of (7) is neither an open-case
certificate nor a Keller counterexample. Other no-finite-pole
binomials outside the retained normal-rescaling and parameter
certificates remain unproved by this note.

## 5. Evidence and remaining frontier

The [checker](verify_central_quotients.py) and
[receipt](receipt-central-quotients.json) provide exact finite
controls for these written all-degree arguments, not an exhaustive
degree or group search. The [source audit](sources/central-quotient-import.json)
records the newly used full faithfulness of Milne Theorem 3.4,
with unchanged retained PDF bytes. The earlier forward-comparison
imports, proofs and receipts are preserved. The
[validation](validation-central-quotients.json) separates scoped
checks from unrelated whole-repository failures.

General nonproduct residual branch arrangements, nonnilpotent
geometric fibre images, parametrizable critical curves outside the
retained certificates, uncovered output geometries, both plane
endpoints and unrestricted HC4/JC2 remain open. Formal verification,
independent replay, external review and novelty priority are not claimed.
