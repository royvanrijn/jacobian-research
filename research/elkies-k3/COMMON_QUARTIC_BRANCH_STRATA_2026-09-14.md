# Branch strata of the fixed-parent common-quartic system

On direct11952 alternate Q80, a genus-one common-quartic solution cannot
have the two sections meet the **same nonzero 2-torsion point at all four
branch values**. This excludes an entire nonconstant-abscissa stratum,
without a coefficient-height bound or an ordinary-node assumption.

This component excludes only the same-branch genus-one stratum. The later
[coefficient collision proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
completes the Q80 polynomial genus-zero exclusions. After the separate
[k2 orientation](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md) and
[k1 cubic reciprocity](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md) proofs,
five genus-one allocations remain open. The fixed-MW17 construction target and
a genus-above-one bound for the full construction remain **OPEN**. In
particular, the genus-zero repeated-quartic boundary is not closed here.
The [previous singularity theorem](COMMON_QUARTIC_SINGULARITY_LOCUS_2026-09-14.md)
retains its independent constant-abscissa exclusion and moving-parent controls.

## 1. A complete factor presentation

Use the fixed short equation `y^2=x^3+A(t)x+B(t)` and write

```
f(x_i)=D*r_i^2,   deg x_i, deg r_i <=4,   deg D=4,
h=x_1-x_2 !=0,   gcd(D,4*A^3+27*B^2)=1.                (1)
```

For now choose a coordinate in which both D and h have degree four.
Put `g=gcd(D,h)`, normalized to be monic, and write

```
D=c*lambda*g*d,       h=c*g*u*v,
r_1=(u*a+v*b)/2,      r_2=(v*b-u*a)/2,
x_2=z,               x_1=z+h.                         (2)
```

Here `g,d,u,v` are monic and `c,lambda` are nonzero rational constants.
The degree strata are exactly

```
k=deg g in {0,...,4},       j=deg u in {0,...,4-k},
deg d=4-k,                 deg v=4-k-j,
deg a<=4-j,                deg b<=k+j,    deg z<=4.     (3)
```

There are fifteen `(k,j)` cases. Their equations for the parent are

```
A=lambda*d*a*b - 3*z^2 - 3*z*h - h^2,
B=D*(v*b-u*a)^2/4 + 2*z^3 + 3*z^2*h + z*h^2
  -lambda*z*d*a*b.                                    (4)
```

All degrees in (4) automatically satisfy `deg A<=8`, `deg B<=12`.
For a fixed parent these are **22 coefficient equations**, not definitions
of replacement parents. There are21 coefficient parameters before the
one-dimensional rescaling of `a,b,lambda`. This count is not an emptiness
argument.

The required open conditions include

```
c*lambda !=0,     gcd(d,u*v)=1,     gcd(v,a)=1,
r_1 !=0,         r_2 !=0,          gcd(D,Delta)=1.       (5)
```

For genus one add that `g*d` is squarefree. For genus zero replace that
condition by `g*d=e^2*q`, with e monic linear and q monic squarefree
quadratic. The possibility `gcd(e,q)!=1` is retained: it gives the triple
root pattern. The literal cover scalar is `c*lambda`; it is not removed
by projectivizing a branch vector. After extracting e the base is
`w^2=c*lambda*q`. A rational point on this conic still has to be supplied.

**Completeness of the presentation.** Subtract (1) to get

```
h*(x_1^2+x_1*x_2+x_2^2+A)=D*(r_1-r_2)*(r_1+r_2).
```

Divide D and h by their monic gcd. Their remaining factors are coprime,
so `h/g` divides `(r_1-r_2)*(r_1+r_2)`. After removing its leading
coefficient c, take `u=gcd(h/c/g,r_1-r_2)` monic and put
`v=h/c/g/u`. Then u divides `r_1-r_2` and v divides `r_1+r_2`, even
when u and v share factors. The definitions of a and b give (2), and
subtraction gives (4). The canonical gcd choice also gives `gcd(v,a)=1`.
Conversely (2), (4) give both identities (1) by direct expansion; (5)
ensures that the prescribed gcd and stratum are the intended ones.
No assumption that h or the node polynomials have simple roots is used.

**Infinity requires only nine predetermined patches.** Fix nine distinct
rational numbers s where the parent fibre is smooth. The checker verifies
`s=1,...,9` for each retained parent, using their existing prime1009 models.
The binary degree-eight form `D*h` cannot vanish at all nine. For one s,
the change `t=s+1/v` and the weighted transforms

```
x_new=v^4*x(s+1/v),   r_new=v^4*r(s+1/v),
D_new=v^4*D(s+1/v),
A_new=v^8*A(s+1/v),   B_new=v^12*B(s+1/v)
```

make D and h have degree four and give smooth infinity. Thus fifteen
factor strata in nine patches cover every admissible rational solution,
including originally constant h, lower-degree abscissas and infinity
contacts. The same proof works for repeated D. This is a finite chart
cover, not an exclusion by testing nine parameter values.

## 2. The factors determine the height matrix

Suppose D is squarefree. On `C:w^2=D`, let `P_i=(x_i,r_i*w)`.
All fibres of the pulled-back surface are irreducible and its holomorphic
Euler characteristic is4. Both points miss zero, so their heights are8.
For the canonical factorization above, put `I=P_1.P_2`. Then

```
I=k+2*j,
<P_1,P_2>=4-k-2*j,
height(P_1-P_2)=8+2*k+4*j,
height(P_1+P_2)=24-2*k-4*j.                            (6)
```

At a common branch root, the parent point is smooth nonzero 2-torsion.
If h has order m, the difference identity makes the orders of
`r_1-r_2` and `r_1+r_2` sum to `m-1`. In the local coordinate w, their
intersection orders are twice those orders plus one. This contributes
one for the root of g and twice its multiplicity in u. Away from the
branch divisor, every intersection occurs over both points of C and
contributes twice its multiplicity in u. The same argument in minimal
infinity coordinates covers that fibre. Points cannot meet the node of
an I1 fibre: locally its smooth surface is `xy=t`, whereas a section
through the node would give order at least two on the left.

Thus I has the stated value. Also
`P_1.P_2+P_1.(-P_2)=2*deg(h)=8`. The section height/intersection formula
on a surface with `chi=4` proves the other three identities. In particular
the possible off-diagonal heights are

| Common branch points k | Possible off-diagonal heights |
|---:|:---|
|0|4, 2, 0, -2, -4|
|1|3, 1, -1, -3|
|2|2, 0, -2|
|3|1, -1|
|4|0|

The determinant is always at least48. These formulas describe actual
solutions if present; the fifteen rows do not assert that solutions exist.
For repeated D, use the general height8 argument from the previous note;
do not reuse `k=deg gcd(D,h)` as a count of branch points.

## 3. The unramified arithmetic input on Q80

The following consequence of retained evidence is needed before halving
a sum of new sections. It cannot be inferred from equality of heights.
For the published determinant948 K3, its good131 reduction has

```
Pic(X_131)=Pic(X)=Z^19 with determinant948,
R_131(T)=(T-131)^19*(T+131)*(T^2+212*T+131^2),
#Br(X_131)=1.                                         (7)
```

These are the exact arithmetic lattice, two-saturation and Artin--Tate
results in [the retained131 proof](../elliptic-curves/rank-jump/THE_LAST_GLOBAL_CLASS_IS_ABSENT.md).
Its independent replay and source-bound certificates are reused, not a
new point count or Selmer computation.

Apply the Kummer/smooth-proper-base-change argument from
[the determinant1092 Brauer theorem, section4](../elliptic-curves/notes/DET1092_SURFACE_BRAUER_TRIVIALITY_2026-09-09.md#4-why-two-good-reductions-rule-out-global-nonconstant-classes)
with `ell=2` and `p=131`. Finite-field Kummer and Hochschild--Serre give

```
H^2(Xbar_131,mu_2)^Frob = H^2(X_131,mu_2)
                       = Pic(X_131)/2.
```

Every geometric restriction of a global Kummer lift is therefore a
global divisor class modulo2. Its Brauer image is algebraic. The full
geometric Picard group is rational, so `H^1(Q,Pic(Xbar))=0`; normalization
at zero gives `Br(X)[2]/Br(Q)[2]=0`. Only2-torsion is asserted; the
residue-characteristic131 part is not treated by this argument.

The [direct11952 construction](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md)
is a fibration on this same Q-surface. The Brauer conclusion is consequently
valid on alternate Q80. The sheaf/Leray proof of
[the unramified Kummer theorem](../elliptic-curves/notes/DET1092_UNRAMIFIED_KUMMER_OBSTRUCTION_2026-09-09.md)
uses precisely this2-primary conclusion, full rational Picard generation
and24I1 fibres. Applying that proof to Q80 gives

```
{norm-square classes in Q(t,theta)^*/squares with all valuations even}
    = delta(E(Q(t))),       theta^3+A*theta+B=0.        (8)
```

All spectral places, including the nodal fibres and infinity, are required.
This deduction does not identify geometric Brauer invariants with
geometric Kummer lifts, and does not apply the Picard-rank-two shortcut.

## 4. Exclusion of the entire k=4 genus-one stratum on Q80

Suppose a solution has k=4. Then `h=c*D` for a nonzero rational c
(renaming the leading scalar). At every branch point the sections select
the same nonzero 2-torsion point. Equation (6) gives

```
Gram(P_1,P_2)=diag(8,8),
U_plus=P_1+P_2,    U_minus=P_1-P_2,
height(U_plus)=height(U_minus)=16.                    (9)
```

For the twist Kummer classes use `alpha_i=D*(x_i-theta)`. At a good
branch place the valuation-parity triple is even at the selected root
and odd at the other two. Identical root choices therefore cancel in
`alpha_1*alpha_2`, whose squareclass is
`beta=(x_1-theta)*(x_2-theta)`. Its norm is `D^2*r_1^2*r_2^2`.
Away from branch places, both local points are defined over an unramified
quadratic extension. The local good/I1 Kummer argument gives even spectral
valuations there too. At infinity the coordinate scaling is a fourth
power. Thus beta satisfies every hypothesis of (8).

There is a generic section T with `beta=delta(T)`. On `L=Q(C)`, the
ordinary Kummer classes of both U_plus and U_minus equal the restriction
of `delta(T)`. Hence there are actual L-rational points

```
R_plus=(U_plus+T)/2,       R_minus=(U_minus+T)/2.        (10)
```

The notation means divisibility in E(L), supplied by Kummer exactness;
it is not formal division of a Mordell--Weil word. Torsion on the pullback
is trivial by the height formula. Conjugating (10) therefore gives
`R_plus+sigma(R_plus)=R_minus+sigma(R_minus)=T` exactly.

Let B_plus and B_minus be their bisection images on the K3. They are
genuine bisections because their anti-traces U_plus, U_minus are nonzero.
For any such image B with point R and anti-trace U,

```
2*B^2=(R+sigma(R))^2=-8+2*(R.sigma(R)),
height(U)=8+2*(R.sigma(R)),
B^2=(height(U)-16)/2.                                 (11)
```

In the first line `R+sigma(R)` denotes the sum of divisors, not elliptic
addition. Equations (9), (11) give `B_plus^2=B_minus^2=0`, so both images
have arithmetic genus one. Their normalization is the given genus-one
curve C, hence they are smooth. Moreover their classes modulo inherited
MW are `(P_1+P_2)/2` and `(P_1-P_2)/2`, which are independent.

This contradicts [the complete smooth-genus-one Q80 exclusion](Q80_ALL_SMOOTH_GENUS_ONE_BISECTIONS_2026-09-13.md).
Indeed the two images have the same trace T and the same square, hence
the same divisor class; even within-pencil injectivity suffices. This
proves the claimed exclusion. It does not assume that the four original
image singularities are nodes or that C has a rational point.

## 5. What remains, and why it is a different arithmetic problem

At any branch point outside g the residues of x_1 and x_2 are distinct
roots of the parent cubic. The third root is `-x_1-x_2`. Thus the cubic
**splits completely over that branch residue field**. In particular,

```
-Delta = [h*(2*x_1+x_2)*(x_1+2*x_2)]^2 mod d.          (12)
```

Let Sigma be the smooth splitting curve of the cubic, with an ordered
pair of distinct roots. Its degree over the t-line is6; the24 I1
inertia transpositions each act as three transpositions on its six sheets.
Riemann--Hurwitz gives `2g(Sigma)-2=-12+24*3`, hence **genus31**.
The branch divisor outside g defines a Q-rational effective divisor of
degree `4-k` on Sigma, with the same residue fields as its t-coordinates.

The k=4 gate excludes degree zero. The later
[reciprocity theorem](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md)
also excludes k=3 on Q80, so a remaining genus-one solution requires such
a divisor of degree2,3 or4. If D is irreducible, the only remaining case
is k=0: the entire quartic branch field must split the cubic. This is a
necessary condition, not an assertion that Sigma has no low-degree points.
Its genus31 does not imply absence of degree-four closed points. Nor does
it bound the genus of a base curve whose branch points merely lie there.

An index-three obstruction on the cubic root curve would exclude a
degree-two or degree-four branch divisor even for one polynomial section.
The preliminary small-prime reductions did not establish that hypothesis:
coordinate changes recover local points in the apparent obstructed charts.
No global index is inferred. Enlarging that probe without a new local or
global argument is not a completed exclusion strategy.

The subsequent [local splitting gate](COMMON_QUARTIC_LOCAL_SPLITTING_GATE_2026-09-14.md)
now supplies actual Q_p points on the full splitting curve for all four
parents at every prime from2 through43. All56 fibres have three distinct
local roots, certified by168 disjoint Hensel balls. Thus this entire prime
panel cannot exclude the needed divisors by an empty local splitting curve.
The global index, low-degree rational divisors and full coefficient incidence
remain unknown; local solubility does not close those gates.

The [checker](scripts/verify_common_quartic_branch_strata.sage) verifies
the universal polynomial identities, all15 degree and height cases,
the nine-patch coverage hypotheses on all four parents, both old controls
through nonleading-coordinate patches, and the arithmetic in (7).
The retained131 and Q80 exclusion certificates are hash-bound proof inputs;
their point counts and billions of comparisons are not repeated. Written
local intersection, cohomological and descent arguments remain mathematical
proofs, not formal proof-assistant verification.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 sage -python research/elkies-k3/scripts/verify_common_quartic_branch_strata.sage
```

The [frozen packet and receipt](../artifacts/generated-results/elkies-k3-common-quartic-branch-strata-v1/)
record a40-CPU-second/4GiB cap. No full coefficient elimination, genus31
point search, new parent rank calculation or positive MW17 solution is
claimed. The next unresolved equation problem is (3)--(5) on the fixed
parents. On Q80 only the five genus-one k=0 allocations remain in
this chart, after the subsequent reciprocity, coefficient collision,
nodal orientation and cubic reciprocity proofs. Repeated D on other parents
remains part of the target.
